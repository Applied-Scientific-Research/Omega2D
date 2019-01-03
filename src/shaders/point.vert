R"(
#version 150

uniform mat4 Projection;
in vec4 position;

void main() {
  // pass through projected coordinates
  gl_PointSize = 2.f;
  gl_Position = Projection * vec4(position.xy, 0.f, 1.f);
}
)"
