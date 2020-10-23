R"(
#version 430
layout(location = 0) uniform mat4 Model;
layout(location = 1) uniform mat4 View;
layout(location = 2) uniform mat4 Projection;
layout (points) in;
layout (triangle_strip, max_vertices = 4) out;
in vec4 colorrad[];
out vec2 txcoord;
out vec4 pcolor;
void main() {
   float rad = colorrad[0].w;
   vec4 pos = (View*Model)*gl_in[0].gl_Position;

   vec4 ppos = Projection*pos;
   // minimum radius
   float fudge = 0.005 * ppos.z;
   float newrad = max(rad, fudge);
   //pcolor = vec4(colorrad[0].x, colorrad[0].y, colorrad[0].z, sqrt(rad / newrad));
   pcolor = vec4(colorrad[0].x, colorrad[0].y, colorrad[0].z, 0.01);

   txcoord = vec2(-1,-1);
   gl_Position = Projection*(pos+newrad*vec4(txcoord,0,0));
   EmitVertex();
   txcoord = vec2( 1,-1);
   gl_Position = Projection*(pos+newrad*vec4(txcoord,0,0));
   EmitVertex();
   txcoord = vec2(-1, 1);
   gl_Position = Projection*(pos+newrad*vec4(txcoord,0,0));
   EmitVertex();
   txcoord = vec2( 1, 1);
   gl_Position = Projection*(pos+newrad*vec4(txcoord,0,0));
   EmitVertex();

   EndPrimitive();
}
)"
