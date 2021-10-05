R"(
#version 150
//
// elempack.vert - set color for an element
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

uniform mat4 Projection;
uniform vec4 pos_color;
uniform vec4 neg_color;
uniform vec4 def_color;
//uniform vec4 back_color;
uniform bool use_def;
//uniform bool use_back;
in vec2 pos;
in float str;
out vec4 base_color;

void main() {
  if (use_def) {
    base_color = def_color;
  } else {
  // flip between base colors based on magnitude of strength
    base_color = abs(str)*(pos_color*step(0.0f, str) + neg_color*step(0.0f, -str));
  }
  // To hide vertices
  //base_color = base_color*(1-use_back) + back_color*use_back;

  // pass 1 vert as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(pos.x, pos.y, 0.f, 1.f);
}
)"
