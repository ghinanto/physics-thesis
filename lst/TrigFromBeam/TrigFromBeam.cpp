/*  This plugin triggers only if the signal comes from the beam and not from
 * other sources (e.g. noise, cosmic rays). L0 trigger just reports and builds
 * hits. L1 trigger builds events. TimeSlice is TriggeredEvent.
 *
 */
#include "TrigFromBeam.h"
#include "Event_Classes.h"
#include "tridas_dataformat.hpp"// needed for fine_time
#include "Geometry.h"
#include "PMTHit.h"
//#include "f_dataformat_clas12.hpp" //needed for fine_time
#include "log.hpp"
#include <algorithm>
#include <cmath>

namespace tridas {
namespace tcpu {

// Returns the row of the 3x3 matrix to which the PMT belongs. Counting starts
// at 0. the indexes of the first row (row == 0) are in the range [0, columns[
// the indexes of the second row (row == 1) are in the range [columns,
// 2*columns[ etc... columns should be equal to geometry.pmts_per_floor rows
// should be equal to geometry.floors_per_tower
int row(PMTHit const& seed, const int columns)
{
  const int seed_PMTID{seed.getAbsPMTID()};

  int row{0};
  while (seed_PMTID < row * columns || seed_PMTID >= (row + 1) * columns) {
    ++row;
  }

  return row;
}

// Returns the column of the 3x3 matrix to which the PMT belongs. Counting
// starts at 0. the indexes of the first column (col == 0) are divisible by the
// number of columns. the indexes of the second column (col == 1) minus 1 are
// divisible by the number of columns. etc... columns should be equal to
// geometry.pmts_per_floor
int col(PMTHit const& seed, const int columns)
{
  const int seed_PMTID{seed.getAbsPMTID()};

  int col{0};
  while ((seed_PMTID - col) % columns) {
    ++col;
  }

  return col;
}

// Two seeds are space correlated if the trace they draw lies at an angle of at
// most 45 degree to the vertical Also because cosmic rays comes from above,
// space correlation requires coherent space AND time order e.g. a seed in row 0
// must be earlier than a seed in row 1
bool spaceCorrelation(PMTHit const& seed, PMTHit const& other_seed,
                      Geometry const& geom)
{
  const int rows{3}; // geom.floors_per_tower;
  const int cols{3}; // geom.pmts_per_floor;
  const int seed_PMTID{seed.getAbsPMTID()};
  assert(0 <= seed_PMTID && seed_PMTID < cols * rows && "PMTID not possible");

  const int seed_row{row(seed, cols)};
  const int other_seed_row{row(other_seed, cols)};
  const int seed_col{col(seed, cols)};
  const int other_seed_col{col(other_seed, cols)};

  // If seed is from a PMT lower (higher) than other_seed, the row index is
  // bigger (smaller) and seed must be older (earlier) Also if the sum and
  // subtraction of row and col index of other_seed are less (greater) or equal
  // to the ones of seed then the two seeds draw an angle to the vertical
  // smaller than 45 degree.
  if (seed_row > other_seed_row
      && seed.IsOlderThan(other_seed)) // space and time order respected
    return other_seed_row + other_seed_col <= seed_row + seed_col
        && other_seed_row - other_seed_col <= seed_row - seed_col;

  if (seed_row < other_seed_row
      && other_seed.IsOlderThan(seed)) // space and time order respected
    return other_seed_row + other_seed_col >= seed_row + seed_col
        && other_seed_row - other_seed_col >= seed_row - seed_col;

  return false; // here if other_seed_row == seed_row
}

// Two seeds are time correlated if the time between them is equal to the time
// spent by the cosmic ray to cover the distance between to PMTs
bool timeCorrelation(PMTHit const& seed, PMTHit const& other_seed,
                     Geometry const& geom)
{
  assert(geom.positions.size() >= 9
         && "Not correct number of PMTs to evaluate time correlation of seeds");
  // The side of the square face of each of the 9 scintillators arranged in the
  // 3x3 matrix i.e. the distance between a PMT and the nearest one in orizontal
  // and vertical directions
  const double front_cell_side{
      std::abs(geom.positions[0].x - geom.positions[1].x)};
  // The diagonal of the square face of each of the 9 scintillators arranged
  // i.e. the distance between a PMT and the nearest one in diagonal direction
  const double front_cell_diagonal{std::sqrt(2) * front_cell_side};
  const double light_speed{2.99792458e8};           // m/s
  const double low_cosmic_speed{0.9 * light_speed}; // m/s

  // Cosmic rays can come from right above or diagonal and in range of speeds of
  // about [0.9c, c] So the shortest time between two seeds is given by a cosmic
  // ray that comes from right above at the speed of light The shortest path
  // possible between two seeds is the side of a single cell The longest time
  // between two seeds is given by a comsic ray coming diagonal at the lowest
  // speed possible. In a 3x3 matrix the longest path possible is going through
  // two cells in diagonal directions
  fine_time shortest_correlation_time_frame{
      static_cast<int>(std::round(front_cell_side / light_speed))};
  fine_time longest_correlation_time_frame{
      static_cast<int>(std::round(2 * front_cell_diagonal / low_cosmic_speed))};
  // Absolute value of the subtraction
  fine_time real_time_frame{
      seed.IsOlderThan(other_seed)
          ? seed.get_fine_time() - other_seed.get_fine_time()
          : other_seed.get_fine_time() - seed.get_fine_time()};

  return shortest_correlation_time_frame < real_time_frame
      && real_time_frame < longest_correlation_time_frame;
}

// Returns the energy of the particle that generated a seed.
// The hits +/- 8 ns from seed belongs to the same particle
// TODO: make half_time_frame a function parameter
double getParticleEnergy(PMTHit const& seed, TriggeredEvent const& tev)
{
  fine_time seed_time{seed.get_fine_time()};
  fine_time half_time_frame{8};

  double particle_energy{seed.getCaliCharge()};

  PMTHit const* this_hit = seed.previous();
  // Go back in time from seed to find the previous hits that fall in a 8 ns
  // time frame
  for (PMTHit const* start_hit = reinterpret_cast<PMTHit const*>(tev.sw_hit_);
       this_hit != start_hit
       && seed_time - this_hit->get_fine_time() <= half_time_frame;
       this_hit = this_hit->previous()) {
    particle_energy += this_hit->getCaliCharge();
  }

  this_hit = seed.next();
  // Go forward in time from seed to find the next hits that fall in a 8 ns time
  // frame
  for (PMTHit const* end_hit = reinterpret_cast<PMTHit const*>(tev.ew_hit_);
       this_hit != end_hit
       && this_hit->get_fine_time() - seed_time <= half_time_frame;
       this_hit = this_hit->next()) {
    particle_energy += this_hit->getCaliCharge();
  }

  return particle_energy;
}

void TrigFromBeam(PluginArgs const& args)
{
  const int id         = args.id;
  EventCollector& evc  = *args.evc;
  Geometry const& geom = *args.geom;

  // Initialise the logger
  tridas::Log::init("TrigFromBeam", tridas::Log::INFO);


  unsigned plug_events{0};

  // for each event
  for (int i{0}; i < evc.used_trig_events(); ++i) {
    TriggeredEvent& tev = *(evc.trig_event(i));
    assert(tev.nseeds_[L1TOTAL_ID] && "Triggered event with no seeds");

    TRIDAS_LOG(INFO) << "In triggered_event " << tev.EventID_ << " there are " << tev.nHit_ << " hits.";

    // search for seeds, i.e. PMTHits of the TriggeredEvent that triggered the
    // L1
    std::vector<PMTHit const*> seeds;

    PMTHit const* end_hit = reinterpret_cast<PMTHit const*>(tev.ew_hit_);
    for (PMTHit const* this_hit = tev.L1_first_seed_; this_hit != end_hit;
         this_hit               = this_hit->next()) {
      assert(this_hit && "invalid current hit address");
      assert(end_hit && "invalid end hit address");

      TRIDAS_LOG(INFO) << "Hit ";
      if (this_hit->isSeed()) {
        seeds.push_back(this_hit);
        TRIDAS_LOG(INFO) << "[SEED] ";
      }
      TRIDAS_LOG(INFO) << *this_hit;
    }
    assert(!seeds.empty() && "Seeds not found");


    // Check if this_seed is from a cosmic ray, i.e. if there are other seeds
    // (~same energy) spacetime correlated to this_seed in other PMTs.
    for (auto const& this_seed : seeds) {
      std::vector<PMTHit const*> cosmic_trace;
      cosmic_trace.reserve(seeds.size());

      std::copy_if(seeds.begin(), seeds.end(), std::back_inserter(cosmic_trace),
                   [&geom, &this_seed](auto const& other_seed) {
                     return this_seed != other_seed
                         && spaceCorrelation(*this_seed, *other_seed, geom)
                         && timeCorrelation(*this_seed, *other_seed, geom);
                   });

      /* Check that are all correlated together in a coherent way to avoid
       * accidental correlations e.g. if the first two are from the same column
       * then the third must be from the same column too. However with only 3
       * columns it would be too difficult to say which is the real one coming
       * from the cosmic ray.
       */

      if (!cosmic_trace.empty()) { // seed is from cosmic ray
        // Add the original seed to fill the trace
        cosmic_trace.push_back(this_seed);

        double cosmic_energy{0.};
        TRIDAS_LOG(INFO) << "Detected cosmic ray:";
        for (auto const& cosmic_seed : cosmic_trace) {
          TRIDAS_LOG(INFO)
              << "\nin PMT[" << cosmic_seed->getAbsPMTID() << "] at time "
              << cosmic_seed->get_fine_time() << " with energy "
              << cosmic_seed->getCaliCharge() << '\n';
          cosmic_energy += cosmic_seed->getCaliCharge();
        }
        TRIDAS_LOG(INFO) << "\nTotal energy is " << cosmic_energy;


      } else { // seed is from beam

        TRIDAS_LOG(INFO) << "Detected particle in PMT["
                          << this_seed->getAbsPMTID() << "] at time "
                          << this_seed->get_fine_time() << " with energy "
                          << getParticleEnergy(*this_seed, tev) << '\n';

        // increment seeds counter of triggered event
        ++tev.plugin_nseeds_[id];
        tev.plugin_ok_ = true;
      }
    }

    if (tev.plugin_ok_) {
      ++plug_events;
    }
  }

  // Set plugin results
  evc.set_stats_for_plugin(id, plug_events);
  TRIDAS_LOG(INFO) << "TrigFromBeam triggered " << evc.stats_for_plugin(id)
                    << " of " << evc.used_trig_events() << " events in TTS "
                    << evc.ts_id();
}

} // namespace tcpu
} // namespace tridas
