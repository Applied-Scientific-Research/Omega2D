R"(
#version 150

uniform float str_scale;
in vec4 base_color;
in vec2 txcoord;
in float strength;
out vec4 frag_color;

void main() {
  // use the core function to draw a fuzzy blob
  float rs = dot(txcoord, txcoord);
  //float s = 1./(1.+16.*rs*rs*rs) - 0.06;
  float s = 1./(1.+16.*rs*rs) - 0.06;

  // scale the color intensity by the particle strength
  frag_color = (str_scale * strength * s) * base_color;
}
)"
