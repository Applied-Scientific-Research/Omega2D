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

uniform float str_scale;
in vec4 base_color;
in vec2 txcoord;
in float strength;
out vec4 frag_color;

void main() {
  // use *some* core function to draw a fuzzy blob
  float rs = dot(txcoord, txcoord);

  // looks more like the compact Gaussian
  float s = 1.0f/(1.0f+16.0f*rs*rs*rs) - 0.06f;

  // looks more like a Gaussian
  //float s = 1.0f/(1.0f+16.0f*rs*rs) - 0.06f;

  // scale the color intensity by the particle strength
  frag_color = (str_scale * strength * s) * base_color;
}
)"
