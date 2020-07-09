R"(
#version 150

in vec4 base_color;
out vec4 frag_color;

void main() {
  // fragment is simply base color
  frag_color = base_color;
  //frag_color = vec4(1.0f, 1.0f, 1.0f, 0.5f);
}
)"
