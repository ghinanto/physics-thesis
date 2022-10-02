/* PT2TXT
 * Take a post-trigger file and extract data to column arranged text.
 *
 */
#include "tridas_dataformat.hpp"
#include <tools/ptfile_reader.hpp>
#include <boost/program_options.hpp>
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>

template<class... Args>
void print(std::ostream& out, Args const&... out_args)
{
  ((out << out_args << " "), ...);
}

template<class KeyType, class ValueType>
void print_keys(std::ostream& out, std::map<KeyType, ValueType> map,
                std::string sep = " ")
{
  std::for_each(map.begin(), map.end(), [&out, &sep](auto const& attribute) {
    KeyType key = attribute.first;
    out << key << sep;
  });
}

template<class KeyType, class ValueType, class UnaryFunction>
void print_keys(std::ostream& out, std::map<KeyType, ValueType> map,
                UnaryFunction on_key, std::string sep = " ")
{
  std::for_each(map.begin(), map.end(),
                [&out, &sep, &on_key](auto const& attribute) {
                  KeyType key = attribute.first;
                  on_key(key);
                  out << key << sep;
                });
}

template<class KeyType, class ValueType>
void print_values(std::ostream& out, std::map<KeyType, ValueType> map,
                  std::string sep = " ")
{
  std::for_each(map.begin(), map.end(), [&out, &sep](auto const& attribute) {
    out << attribute.second << sep;
  });
}

void progress(unsigned long& index, unsigned long const& tot,
              std::ostream& out_interface)
{
  double progress_double =
      static_cast<double>(index++) / static_cast<double>(tot) * 100.;
  unsigned long progress{0};

  // Write to output only if necessary
  if (progress < static_cast<unsigned long>(progress_double)) {
    progress = static_cast<unsigned long>(progress_double);
    out_interface << '\r' << progress << '%';
    out_interface.flush();
  }
}

namespace po = boost::program_options;
typedef tridas::post_trigger::sample::uncompressed Sample;

int main(int argc, char** argv)
{
  std::string inptfile;
  std::string outtxtfile;

  po::options_description desc("Options");
  // clang-format off
  desc.add_options()
  ("help,h", "Show this message.")
  ("input,i", po::value<std::string>(&inptfile)->required(), 
  "The post-trigger file to read.")
  ("output,o", po::value<std::string>(&outtxtfile)->required(), 
  "The text file to write.")
  ;
  // clang-format on

  try {
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return EXIT_SUCCESS;
    }

    po::notify(vm);
  } catch (po::error const& e) {
    std::cerr << "[PT2TXT] Error: " << e.what() << '\n';
    std::cerr << desc << std::endl;
    return EXIT_FAILURE;
  } catch (std::exception const& e) {
    std::cerr << "[PT2TXT] Error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }

  // std::ostream& out = std::cout;
  std::ofstream outtxtfile_stream(outtxtfile);
  std::ostream& out           = outtxtfile_stream;
  std::ostream& out_interface = std::cout;
  bool first_line{true};

  // extract datacard from pt file and check plugins

  tridas::post_trigger::PtFileReader<Sample> const ptreader(inptfile);
  std::map<std::string, int> pt_attributes = {
      {"run_number", ptreader.runNumber()},
      {"file_number", ptreader.fileNumber()},
      {"run_start_time", ptreader.runStartTime()},
      {"datacard_size", ptreader.datacardSize()},
      {"number_of_events", ptreader.nEvents()},
      {"number_of_timeslices", ptreader.nTS()},
      {"effective_file_size", ptreader.fileSize()}};

  unsigned long nTS = ptreader.nTS();
  // out_interface << "TimeSlice[" << nTS << "]:\n";
  unsigned long ts_count{0};
  for (auto const& ts : ptreader) {
    std::map<std::string, int> ts_attributes = {
        {"TS_ID", ts.id()}, {"NEvents", ts.nEvents()}, {"TS_size", ts.size()}};

    unsigned long te_count{0};
    for (auto const& te : ts) {
      std::map<std::string, int> te_attributes = {
          {"EventTag", te.eventTag()},
          {"EventID", te.id()},
          {"EventL", te.size()},
          {"nHit", te.nHit()},
          {"StartTime", te.startTime().count()},
          {"TSCompleted", te.tsComplete()},
          {"nseeds[0]", te.nseeds(0)},
          {"plugin_trigtype[0]", te.plugin_trigtype(0)},
          {"plugin_nseeds[0]", te.plugin_nseeds(0)}};

      unsigned long hit_count{0};
      for (auto const& hit : te) {
        std::map<std::string, int> hit_attributes = {
            {"nFrames", hit.nFrames()}};

        for (long i{0}; i < hit.nFrames(); ++i) {
          auto& dfheader = hit.frameHeader(i++);

          std::map<std::string, int> df_attributes = {
              {"Charge", dfheader.Charge},
              {"EFCMID", dfheader.EFCMID},
              {"FrameCounter", dfheader.FrameCounter},
              {"PMTID", dfheader.PMTID},
              {"T1ns", dfheader.T1ns},
              {"TowerID", dfheader.TowerID}};

          if (first_line) {
            print_keys(out, df_attributes);
            print_keys(out, hit_attributes);
            print_keys(out, te_attributes, [](std::string& key) {
              if (key == "StartTime") {
                key +=
                    ("[" + std::to_string(decltype(te.startTime())::period::num)
                     + "/"
                     + std::to_string(decltype(te.startTime())::period::den)
                     + "]sec");
              }
            });
            print_keys(out, ts_attributes);
            print_keys(out, pt_attributes);
            out << '\n';
            first_line = false;
          }

          print_values(out, df_attributes);
          print_values(out, hit_attributes);
          print_values(out, te_attributes);
          print_values(out, ts_attributes);
          print_values(out, pt_attributes);
          out << '\n';

          // std::cout << count
        }
        // std::cout << "\nHit[" << te.nHit() << "]:" << hit_count++;
      }
      // std::cout.flush();
      // std::cout << "\nEvent[" << ts.nEvents() << "]:" << te_count++;
    }
    // std::cout.flush();
    progress(ts_count, nTS, out_interface);
  }

  return EXIT_SUCCESS;
}