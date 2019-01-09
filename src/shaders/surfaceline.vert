R"(
#version 150

uniform mat4 Projection;
uniform vec4 pos_color;
uniform vec4 neg_color;
uniform vec4 def_color;
in float px;
in float py;
in float rawstr;
out vec4 base_color;
out float strength;

void main() {
  // flip between base colors based on magnitude of strength
  //base_color = pos_color*step(0.0f, rawstr) + neg_color*step(0.0f, -rawstr);
  base_color = def_color;

  // magnitude of strength scales fragments
  //strength = 1.f;
  strength = abs(rawstr);

  // make 4 verts as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(px, py, 0.f, 1.f);
}
)"
