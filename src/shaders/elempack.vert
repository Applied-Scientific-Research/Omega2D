R"(
#version 150

uniform mat4 Projection;
uniform vec4 pos_color;
uniform vec4 neg_color;
uniform vec4 def_color;
uniform bool one_color;
in vec2 pos;
in float str;
out vec4 base_color;

void main() {
  // flip between base colors based on magnitude of strength
  base_color = pos_color*step(0.0f, str) + neg_color*step(0.0f, -str);
  if (one_color) {
    base_color = def_color;
  }
  // scale the default color by the input strength
  //base_color = str * def_color;

  // pass 1 vert as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(pos.x, pos.y, 0.f, 1.f);
}
)"
