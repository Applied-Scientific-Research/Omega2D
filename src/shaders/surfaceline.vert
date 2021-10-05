R"(
#version 150
//
// surfaceline.vert - project a line
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
in float px;
in float py;
in float rawstr;
out vec4 base_color;
out float strength;

void main() {
  // flip between base colors based on magnitude of strength
  //base_color = pos_color*step(0.0f, rawstr) + neg_color*step(0.0f, -rawstr);
  base_color = def_color;

  // magnitude of strength scales fragments
  //strength = 1.f;
  strength = abs(rawstr);

  // make 4 verts as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(px, py, 0.f, 1.f);
}
)"
