R"(
#version 150

uniform mat4 Projection;
uniform vec4 pos_color;
uniform vec4 neg_color;
uniform float rad_scale;
in vec4 quad_attr;
in float px;
in float py;
in float rad;
in float sx;
out vec4 base_color;
out vec2 txcoord;
out float strength;

void main() {
  // flip between base colors based on whether strength points in or out of the screen
  base_color = pos_color*step(0.0f, sx) + neg_color*step(0.0f, -sx);

  // how large is the spot, actually
  float rscale = 2.5f * rad * rad_scale;

  // if particle would be too small, make it dimmer instead
  // use top-left of Proj matrix for scaling
  float draw_scale = max(0.002/Projection[0][0], rscale);
  float scale_ratio = rscale / (rad*draw_scale);

  // pass through texture coordinates
  txcoord = quad_attr.xy;

  // magnitude of strength density scales fragment color
  strength = abs(sx)*scale_ratio*scale_ratio;

  // make 4 verts as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(px + draw_scale*quad_attr.x, py + draw_scale*quad_attr.y, 0.f, 1.f);
}
)"
