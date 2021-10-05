/*
 * JsonHelper.h - Methods to help with reading and writing JSON files
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
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
