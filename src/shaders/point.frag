R"(
#version 150

in vec4 base_color;
in vec2 txcoord;
out vec4 frag_color;

void main() {
  // use *some* core function to draw a fuzzy blob
  float rs = dot(txcoord, txcoord);
  float s = 1./(1.+16.*rs*rs*rs) - 0.06;

  // scale the color intensity by the particle strength
  frag_color = s * base_color;
}
)"
