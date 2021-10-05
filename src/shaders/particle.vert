R"(
#version 150
//
// particle.vert - scale and project a particle
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
uniform float rad_scale;
in vec4 quad_attr;
in float px;
in float py;
in float rad;
in float sx;
out vec4 base_color;
out vec2 txcoord;
out float strength;

void main() {
  // flip between base colors based on whether strength points in or out of the screen
  base_color = pos_color*step(0.0f, sx) + neg_color*step(0.0f, -sx);

  // how large is the spot, actually
  float rscale = 2.5f * rad * rad_scale;

  // if particle would be too small, make it dimmer instead
  // use top-left of Proj matrix for scaling
  float draw_scale = max(0.002/Projection[0][0], rscale);
  float scale_ratio = rscale / (rad*draw_scale);

  // pass through texture coordinates
  txcoord = quad_attr.xy;

  // magnitude of strength density scales fragment color
  strength = abs(sx)*scale_ratio*scale_ratio;

  // make 4 verts as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(px + draw_scale*quad_attr.x, py + draw_scale*quad_attr.y, 0.f, 1.f);
}
)"
