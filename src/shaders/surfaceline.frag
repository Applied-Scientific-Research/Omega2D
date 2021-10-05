R"(
#version 150
//
// surfaceline.frag - color pass-through
//
// (c)2019 Applied Scientific Research, Inc.
//         Mark J. Stock <markjstock@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

uniform float str_scale;
in vec4 base_color;
in float strength;
out vec4 frag_color;

void main() {
  // scale the color intensity by the particle strength
  //frag_color = (str_scale * strength) * base_color;
  frag_color = base_color;
  //frag_color = vec4(1.0f, 1.0f, 1.0f, 0.5f);
}
)"
