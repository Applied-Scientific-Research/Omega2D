/*
 * JsonHelper.h - Methods to help with reading and writing JSON files
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Simulation.h"
#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "MeasureFeature.h"
#include "RenderParams.h"

#include <string>

nlohmann::json read_json(const std::string filename);

void parse_json(Simulation&,
                std::vector<std::unique_ptr<FlowFeature>>&,
                std::vector<std::unique_ptr<BoundaryFeature>>&,
                std::vector<std::unique_ptr<MeasureFeature>>&,
                RenderParams&,
                nlohmann::json);

void write_json(Simulation&,
                std::vector<std::unique_ptr<FlowFeature>>&,
                std::vector<std::unique_ptr<BoundaryFeature>>&,
                std::vector<std::unique_ptr<MeasureFeature>>&,
                const RenderParams&,
                const std::string);
