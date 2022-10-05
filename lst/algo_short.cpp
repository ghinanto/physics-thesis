for (auto const& this_seed : seeds) {
    std::vector<PMTHit const*> cosmic_trace;
    cosmic_trace.reserve(seeds.size());

    std::copy_if(seeds.begin(), seeds.end(), std::back_inserter(cosmic_trace),
           [&geom, &this_seed](auto const& other_seed) {
             return this_seed != other_seed
             && spaceCorrelation(*this_seed, *other_seed, geom)
             && timeCorrelation(*this_seed, *other_seed, geom);
           });

    if (!cosmic_trace.empty()) { // seed is from cosmic ray
        // get information about cosmic ray
    } else { // seed is from beam
        // select seed to be saved in the PT file
    }
}

