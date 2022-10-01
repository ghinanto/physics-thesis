#include "TowerTimeSlice.h"
#include "log.hpp"
#include "pmt_numbering.hpp"
#include "trig_sim_details.hpp"
#include <tools/ptfile_reader.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <iostream>

// For EM capabilities
#include <DAQ/EM/inc/FileWriter.h>

// For TSV capabilities
#include "inter_objects_communication/producer_consumer.hpp"

// For HM configuration
#include "HMParameters.h"

#include "subprocess.hpp"

typedef tridas::post_trigger::sample::uncompressed Sample;

// Copied from RunTSV.cpp
struct DoneTSMessage
{
  TS_t timeslice_id;
  TcpuId tcpu_id;
};

namespace tpt              = tridas::post_trigger;
namespace po               = boost::program_options;
namespace em               = tridas::em;
namespace strict_numbering = tridas::strict_numbering;

using namespace tridas::femsim;

int main(int argc, char* argv[])
{
  int const n_sectors = 1;
  std::string inptfile;
  std::string outptfile;
  std::string tcpu_path;
  std::string cardfile;

  po::options_description desc("Options");

  // clang-format off
  desc.add_options()
    ("help,h", "Print help messages.")
    ("tcpu",
     po::value<std::string>(&tcpu_path)
         ->required()
         ->value_name("filename"),
     "Path of the RunTCPU executable.")
    ("datacard-input,di",
     po::value<std::string>(&cardfile)
         ->required()
         ->value_name("filename"),
     "Path of the datacard to read.")
    ("input,i",
     po::value<std::string>(&inptfile)
         ->required()
         ->value_name("filename"),
     "Input Post Trigger file.")
    ("output,o",
     po::value<std::string>(&outptfile)
         ->required()
         ->value_name("filename"),
     "Output Post Trigger file.");
  // clang-format on

  try {
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return EXIT_SUCCESS;
    }

    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << "TrigSim: Error: " << e.what() << '\n';
    std::cerr << desc << std::endl;
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cerr << "TrigSim: Error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }

  // Initialise the logger
  tridas::Log::init("TrigSim", tridas::Log::INFO);

  // Create the datacard
  Configuration conf;
  read_json(cardfile, conf);

  int const pmts_per_floor   = conf.get<int>("DETECTOR_GEOMETRY.PMTS");
  int const floors_per_tower = conf.get<int>("DETECTOR_GEOMETRY.FLOORS");
  int const towers           = conf.get<int>("DETECTOR_GEOMETRY.TOWERS");
  int const floor_offset     = conf.get<int>("DETECTOR_GEOMETRY.FLOOR_OFFSET");
  int const socket_buffer_size =
      conf.get<int>("INTERNAL_SW_PARAMETERS.OS_SOCKET_BUFFER_SIZE");

  // Prepare the input data
  strict_numbering::PmtId pmt_id(floors_per_tower,
                                 pmts_per_floor); // 5 floors, 16 pmt per floor

  int const n_tcpu_threads = 1;
  // Launch the TCPU
  std::string const tcpu_cmd_line = tcpu_path + " -i 0 -d " + cardfile;
  //+ " 2>/tmp/stderr.txt 1>/tmp/stdout.txt";
  Subprocess tcpu_process(tcpu_cmd_line);

  // Create connections
  network::Context context;
  network::Consumer tsv2tcpu(context, socket_buffer_size);
  std::vector<network::ProducerPtr> hm2tcpu;
  network::Consumer tcpu2em(context, socket_buffer_size);

  tsv2tcpu.tcpConnect("localhost",
                      conf.get<std::string>("TCPU.BASE_CTRL_PORT"));

  for (int i = 0; i < n_sectors; ++i) {
    hm2tcpu.push_back(network::ProducerPtr(
        new network::Producer(context, socket_buffer_size)));
    hm2tcpu[i]->tcpConnect("localhost",
                           conf.get<std::string>("TCPU.BASE_DATA_PORT"));
  }

  tcpu2em.tcpBind(conf.get<std::string>("EM.DATA_PORT"));

  em::FileWriter filewriter(
      tcpu2em, conf.get<int>("INTERNAL_SW_PARAMETERS.DELTA_TS"), "localhost",
      9999, 1, time(0), std::numeric_limits<std::size_t>::max(), cardfile,
      outptfile);

  boost::thread em_th(boost::ref(filewriter));

  tpt::PtFileReader<Sample> const reader(inptfile);

  int const nts       = reader.nTS();
  int const step      = nts < 10 ? 1 : nts / 10;
  int ts_count        = 0;
  int const npmt_sect = towers * floors_per_tower * pmts_per_floor / n_sectors;

  TRIDAS_LOG(INFO) << "Starting!\n";

  std::list<tpt::TimeSlice<Sample>> pt_file;
  for (tpt::PtFileReader<Sample>::const_iterator it = reader.begin(),
                                                 et = reader.end();
       it != et; ++it) {
    pt_file.push_back(*it);
  }
  // Sort timeslices by timeslice id
  pt_file.sort([](tpt::TimeSlice<Sample> const& first,
                  tpt::TimeSlice<Sample> const& second) {
    return first.id() < second.id();
  });

  for (auto const& ts : pt_file) {
    // convert TS to a number of STS
    std::map<int, std::vector<unsigned char>> const tts = convert(ts, pmt_id);

    TRIDAS_LOG(INFO) << "Number of PMTs in next TS: " << tts.size() << '.';

    // The following is just to block until a token has been received.
    {
      DoneTSMessage tsmessage;
      tsv2tcpu.timedRecv(&tsmessage, sizeof(tsmessage), -1);
    }

    TRIDAS_LOG(INFO) << "Sending TS " << ts.id() << '.';

    for (int i = 0; i < n_sectors; ++i) {
      std::map<int, std::vector<unsigned char>>::const_iterator const lower =
          tts.lower_bound(i * npmt_sect);

      std::map<int, std::vector<unsigned char>>::const_iterator const upper =
          tts.upper_bound((i + 1) * npmt_sect);

      if (i == n_sectors - 1)
        assert(upper == tts.end());

      std ::size_t sent_size =
          send(lower, upper, ts.id(), i, npmt_sect, *hm2tcpu[i]);

      TRIDAS_LOG(INFO) << "Sent " << sent_size << " Bytes";
    }

    TRIDAS_LOG(INFO) << "TS " << ts.id() << " sent.";

    ++ts_count;
    if (ts_count % step == 0) {
      double const percentage = (100. * ts_count) / nts;
      TRIDAS_LOG(INFO) << "Work status: " << std::setprecision(2) << percentage;
    }
  }

  // The following is just to block until a token has been received.
  for (int i = 0; i < n_tcpu_threads; ++i) {
    DoneTSMessage tsmessage;
    tsv2tcpu.timedRecv(&tsmessage, sizeof(tsmessage), -1);
  }

  tcpu_process.kill(2);
  filewriter.stop();
  em_th.join();
}
