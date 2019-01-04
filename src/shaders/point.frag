R"(
#version 150

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
