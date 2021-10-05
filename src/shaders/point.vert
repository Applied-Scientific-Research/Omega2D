R"(
#version 150
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
uniform vec4 def_color;
uniform float rad;
in vec4 quad_attr;
in float px;
in float py;
out vec4 base_color;
out vec2 txcoord;

void main() {
  // color pass-through
  base_color = def_color;

  // if particle would be too small, make it dimmer instead
  // use top-left of Proj matrix for scaling
  float draw_rad = max(0.002/Projection[0][0], rad);
  base_color = base_color * pow(rad / draw_rad, 2);

  // pass through texture coordinates
  txcoord = quad_attr.xy;

  // make 4 verts as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(px + draw_rad*quad_attr.x, py + draw_rad*quad_attr.y, 0.f, 1.f);
}
)"
