R"(
#version 150

uniform mat4 Projection;
uniform vec4 pos_color;
uniform vec4 neg_color;
in vec2 position;
in float rawstr;
out vec4 base_color;
out float strength;

void main() {
  // flip between base colors based on magnitude of strength
  base_color = pos_color*step(0.0f, rawstr) + neg_color*step(0.0f, -rawstr);
  //base_color = pos_color;

  // magnitude of strength scales fragments
  //strength = abs(position.z);
  //strength = 1.f;
  strength = abs(rawstr);

  // make 4 verts as a single primitive and set texture coords - see other shaders
  //float rad = position.w;
  //float rad = 0.005f;
  //gl_Position = Projection * vec4(position.xy + 2.5f*rad*quad_attr.xy, 0.f, 1.f);
  //gl_PointSize = 5.f;
  gl_Position = Projection * vec4(position.xy, 0.f, 1.f);
}
)"
