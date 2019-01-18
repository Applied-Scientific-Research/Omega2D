R"(
#version 150

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
