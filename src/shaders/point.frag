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

in vec4 base_color;
in vec2 txcoord;
out vec4 frag_color;

void main() {
  // use *some* core function to draw a fuzzy blob
  float rs = dot(txcoord, txcoord);
  float s = 1.0f/(1.0f+16.0f*rs*rs*rs) - 0.06f;

  // scale the color intensity by the particle strength
  frag_color = s * base_color;
}
)"
