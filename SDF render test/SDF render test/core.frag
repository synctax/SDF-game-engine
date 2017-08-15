#version 330 core

//#define WIDTH (500)
//#define WIDTHINV (0.002)
//#define HEIGHT (500)
//#define HEIGHTINV (0.002)

#define WIDTH (320)
#define WIDTHINV (0.003125)
#define HEIGHT (800)
#define HEIGHTINV (0.00125)

//#define WIDTH (1280)
//#define WIDTHINV (0.00078125)
//#define HEIGHT (800)
//#define HEIGHTINV (0.00125)

#define FOV_X (60)
#define FOV_Y (60)
#define MAX_ITERATIONS (300)
#define T_MAX (1000)
#define T_MIN (0)
#define FORCE_HIT (0)
#define epsilon (0.0001)
#define PI 3.14159265
#define HASHSCALE3 vec3(.1031, .1030, .0973)


///////////////////////////  LIBRARY FUNCTIONS  ///////////////////////////
/*         PRIMITIVES        */
float fSphere( vec3 p, vec3 args );
float fRoundBox( vec3 p, vec3 b);
float fBox( vec3 p, vec3 b );
float fTorus( vec3 p, vec3 args );
float fCylinder( vec3 p, vec3 c );
float fCone( vec3 p, vec3 args );
float fPlane( vec3 p, vec3 n );
float fHexPrism( vec3 p, vec3 args );
float fTriPrism( vec3 p, vec3 args );
float fCapsule(vec3 p, vec3 args);
float fCappedCylinder( vec3 p, vec3 args );
float fCappedCone(vec3 p, vec3 c );
float fEllipsoid( vec3 p, vec3 r );



/*         DOMAIN MANIPULATION OPERATORS        */
void pR(inout vec2 p, vec3 args);
void pR45(inout vec2 p);
float pMod1(inout float p, vec3 args);
float pModPolar(inout vec2 p, vec3 args);
vec2 pMod2(inout vec2 p, vec3 args);
vec3 pMod3(inout vec3 p, vec3 size);



/*         OBJECT COMBINATION OPERATORS        */
vec2 fOpUnion(vec2 a, vec2 b);
vec2 fOpIntersect(vec2 a, vec2 b);
vec2 fOpDifference(vec2 a, vec2 b);
vec2 fOpUnionChamfer(vec2 a, vec2 b, vec3 args);
vec2 fOpIntersectionChamfer(vec2 a, vec2 b, vec3 args);
vec2 fOpDifferenceChamfer (vec2 a, vec2 b, vec3 args);
vec2 fOpUnionRound(vec2 a, vec2 b, vec3 args);
vec2 fOpIntersectionRound(vec2 a, vec2 b, vec3 args);
vec2 fOpDifferenceRound (vec2 a, vec2 b, vec3 args);

/*                     MISC                    */
float noise(vec3 x);
vec2 rand(vec2 p);
vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n );
vec3 texBump( sampler2D tx, in vec3 p, in vec3 n, float bf);

/*         OBJECT MODIFICATION OPERATORS       */
vec2 fBumpify(vec3 p, vec2 shape, vec3 args);

//////////////////////////////////////////////////////////////////////////////////

out vec4 color;
in vec2 TexCoord;
uniform vec3 playerPosition;
uniform vec3 playerView;
uniform sampler2D redSandTex;
uniform sampler2D redSandBMP;

const float pixelWidth = 2*WIDTHINV;
const float pixelHeight = 2*HEIGHTINV;
const float xOffset = -1;
const float yOffset = -1;
const float aspect = WIDTH*HEIGHTINV;
const vec3 ambientLight = vec3(0,0,0);
uniform int sceneTree[50];

struct sdfNode{
    int nodeType;
    int subType;
    vec3 args;
    int extra;
    int flags;
    int children[2];
    int addr;
};

struct pointLight{
    vec3 pos;
    vec3 col;
    float intensity;
};

struct Material{
    vec3 color;
    float reflectivity;
    float specular;
//    bool hasTex;
//    sampler2D tex;
//    bool hasBMP;
};

Material materials[2] = Material[2]( Material(vec3(0.49,0.21176,0.133333),0.0,0.0),
                                    Material(vec3(0.6784,0.549,0.3647),0.1,0.1));

const int numLights = 2;
pointLight lights[numLights] = pointLight[numLights](pointLight(vec3(100,100,50), vec3(1,1,1), 1.0),
                                                    pointLight(vec3(-100,100,-50), vec3(1,1,1), 1.0));

vec2 makePath(vec3 p){
    vec2 path = vec2(fBox(p,vec3(5,0.25,100)),1);
    p +=vec3(55,0,-95);
    path = fOpUnion(path,vec2(fBox(p,vec3(50,0.25,5)),1));
    p += vec3(-50,0,0);
    pR45(p.xz);
    path = fOpUnion(path,vec2(fBox(p+vec3(-80,0,3.25),vec3(75,0.25,5)),1));
    pR45(p.zx);
    path = fOpUnion(path,vec2(fBox(p+vec3(-118.35,0,-105.83),vec3(10,0.25,5)),1));
    return path;
}

vec2 scene(vec3 p){
    vec2 ground = vec2(p.y,0);
    vec2 path = makePath(p);
    return fOpUnion(ground,path);
}

vec3 getSceneNormal( vec3 pos, vec3 ray, float distance)
{
    
    vec2 d = vec2(-1,1) * .5 * distance / WIDTH;
    vec3 p0 = pos+d.xxx;
    vec3 p1 = pos+d.xyy;
    vec3 p2 = pos+d.yxy;
    vec3 p3 = pos+d.yyx;
    
    float f0 = scene(p0).x;
    float f1 = scene(p1).x;
    float f2 = scene(p2).x;
    float f3 = scene(p3).x;
    
    vec3 grad = p0*f0+p1*f1+p2*f2+p3*f3 - pos*(f0+f1+f2+f3);
    
    float gdr = dot ( grad, ray );
    grad -= max(.0,gdr)*ray;
    
    return normalize(grad);
}

int enhancedSphereTrace(out vec3 p, vec3 o, vec3 d){
    float omega = 1.2;
    float t = T_MIN;
    float candidate_error = 1000000000;
    float candidate_t = T_MIN;
    float previousRadius = 0;
    float stepLength = 0;
    float pixelRadius = WIDTHINV;
    int i;
    for (i = 0; i < MAX_ITERATIONS; ++i) {
        float radius = scene(d*t + o).x;
        bool sorFail = omega > 1 && (radius + previousRadius) < stepLength;
        if (sorFail) {
            stepLength -= omega* stepLength;
            omega = 1;
        } else {
            stepLength = radius*omega;
        }
        previousRadius = radius;
        
        float error = radius / t;
        if (!sorFail && error < candidate_error) {
            candidate_t = t;
            candidate_error = error;
        }
        if (!sorFail && error < pixelRadius || t > T_MAX ) break;
        t += stepLength;
    }
    if ((t > T_MAX || candidate_error > pixelRadius) && FORCE_HIT==0) return 0;
    
    p = o + (candidate_t*d);
    return 1;
}

float softShadow( vec3 ro, vec3 rd, float maxt, float k )
{
    float res = 1.0;
    for( float t=0; t < maxt; )
    {
        float h = scene(ro + rd*t).x;
        if( h<0.01 )
            return 0.0;
        res = min( res, k*h/t );
        t += h;
    }
    return res;
}
float calculateAO(in vec3 pos, in vec3 nor, float range)
{
    float sca = 2.0, occ = 0.0;
    for( int i=0; i<5; i++ ){
        
        float hr = 0.01 + float(i)*0.5/4.0;
        float dd = scene(nor * hr + pos).x;
        occ += (hr - dd)*sca;
        sca *= 0.7;
    }
    return clamp( 1.0 - occ/range, 0.0, 1.0 );
}

vec3 sky( vec3 ray )
{
    return mix( vec3(.8), vec3(0), exp2(-(1.0/max(ray.y,.01))*vec3(.4,.6,1.0)) );
}


vec3 shade( vec3 pos, vec3 ray, vec3 normal, vec3 lightDir, pointLight light, Material mat,float distance )
{
    normal = texBump(redSandTex,pos/100,normal,0.1);
    vec3 ambient = mix( vec3(.05,.03,.04), vec3(.1), (-normal.y+1.0) );
    float aoRange = distance/20.0;
    //float occlusion = max( 0.0, 1.0 - scene( pos + normal*aoRange ).x/aoRange );
    float occlusion = calculateAO(pos,normal,aoRange);
    occlusion = exp2( -3.0*pow(occlusion,2.0) );
    ambient *= occlusion;
    
    vec3 l = light.col*max(.0,dot(normal,lightDir));
    l += ambient;
    
    vec3 h = normalize(lightDir-ray);
    vec3 specular = light.col*pow(max(.0,dot(normal,h)), 100.0)*100.0/32.0;
    specular *= mat.specular;
    
    vec3 rray = reflect(ray,normal);
    vec3 p, reflection = vec3(1);
    //    if (enhancedSphereTrace(p,pos,ray) == 1){
    //        reflection = materials[int(scene(p).y)].color;
    //    }else{
    //        reflection = sky(rray);
    //    }
    reflection *= mat.reflectivity;
    reflection *= occlusion;
    occlusion = max( 0.0, 1.0 - scene( pos + rray*aoRange ).x/(aoRange*dot(rray,normal)) );
    occlusion = exp2( -2.0*pow(occlusion,2.0) );
    reflection *= occlusion*mat.reflectivity;
    
    float fresnel = pow( 1.0+dot(normal,ray), 5.0 );
    fresnel = mix( .0, .5, fresnel );
    
    vec3 albedo = mat.color;
    
    return mix( l*albedo, reflection, fresnel ) + specular;
}

vec3 colorAtPoint(vec3 p, vec3 rayDir,float dist){
    Material mat = materials[int(scene(p).y)];
    vec3 R, n = getSceneNormal(p,rayDir,dist), lightDir;
    float dist2;
    pointLight light;
    for (int i = 0; i < numLights; i++){
        light = lights[i];
        lightDir = light.pos-p; dist2 = length(lightDir);
        lightDir = normalize(lightDir);
        R+=shade(p,rayDir,n,lightDir,light,mat,dist)*light.intensity;
    }
    return R;
}

vec3 rayToScreenPoint(vec3 o, vec3 d, float x, float y){
    float xPos = ((pixelWidth*x)+xOffset)*aspect;
    float yPos =(pixelHeight*y)+yOffset;
    
    vec3 xHat = normalize(cross(vec3(0,1,0),d))*1.2;
    vec3 yHat = normalize(cross(xHat, d));
    
    return normalize(xPos*xHat+yPos*yHat+d);
}

vec3 trace(vec3 dir, vec3 pos){
    vec3 col;
    vec3 p;
    if (enhancedSphereTrace(p,pos,dir) == 1){
        col = colorAtPoint(p,dir,distance(p,pos));
    }else{
        col = sky(dir);
        for (int i = 0; i < numLights; i++){
            vec3 lightDir = normalize(lights[i].pos-pos);
        }
        
    }
    return col;
}

void main(){
    float tx = TexCoord.x*WIDTH;
    float ty = TexCoord.y*HEIGHT;
    vec3 dir = rayToScreenPoint(playerPosition, playerView, tx, ty);
    
    color = vec4(trace(dir,playerPosition),1);
    //color = vec4(texture(redSandBMP,TexCoord));
}
























////////////////////////////////////////////////////////////////
//
//                           HG_SDF
//
//     GLSL LIBRARY FOR BUILDING SIGNED DISTANCE BOUNDS
//
//     version 2015-12-15 (initial release)
//
//     Check http://mercury.sexy/hg_sdf for updates
//     and usage examples. Send feedback to spheretracing@mercury.sexy.
//
//     Brought to you by MERCURY http://mercury.sexy
//
//
//
// Released as Creative Commons Attribution-NonCommercial (CC BY-NC)
//
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//
//             HELPER FUNCTIONS/MACROS
//
////////////////////////////////////////////////////////////////

#define PI 3.14159265
#define TAU (2*PI)
#define PHI (1.618033988749895)
// PHI (sqrt(5)*0.5 + 0.5)

// Clamp to [0,1] - this operation is free under certain circumstances.
// For further information see
// http://www.humus.name/Articles/Persson_LowLevelThinking.pdf and
// http://www.humus.name/Articles/Persson_LowlevelShaderOptimization.pdf
#define saturate(x) clamp(x, 0., 1.)

// Sign function that doesn't return 0
float sgn(float x) {
    return (x<0.)?-1.:1.;
}

float square (float x) {
    return x*x;
}

vec2 square (vec2 x) {
    return x*x;
}

vec3 square (vec3 x) {
    return x*x;
}

float lengthSqr(vec3 x) {
    return dot(x, x);
}


// Maximum/minumum elements of a vector
float vmax(vec2 v) {
    return max(v.x, v.y);
}

float vmax(vec3 v) {
    return max(max(v.x, v.y), v.z);
}

float vmax(vec4 v) {
    return max(max(v.x, v.y), max(v.z, v.w));
}

float vmin(vec2 v) {
    return min(v.x, v.y);
}

float vmin(vec3 v) {
    return min(min(v.x, v.y), v.z);
}

float vmin(vec4 v) {
    return min(min(v.x, v.y), min(v.z, v.w));
}




////////////////////////////////////////////////////////////////
//
//             PRIMITIVE DISTANCE FUNCTIONS
//
////////////////////////////////////////////////////////////////
//
// Conventions:
//
// Everything that is a distance function is called fSomething.
// The first argument is always a point in 2 or 3-space called <p>.
// Unless otherwise noted, (if the object has an intrinsic "up"
// side or direction) the y axis is "up" and the object is
// centered at the origin.
//
////////////////////////////////////////////////////////////////

//Sphere: A sphere with radius in args.x
float fSphere( vec3 p, vec3 args ){
    return length(p)-args.x;
}

//Round Box: a rounded box of dimensions b
float fRoundBox( vec3 p, vec3 b){
    return length(max(abs(p)-b,0.0))-0.5;
}

//Box: a box of dimensions b
float fBox( vec3 p, vec3 b ){
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

//Torus: a torus of args.x smallR and args.y bigR
float fTorus( vec3 p, vec3 args ){
    vec2 q = vec2(length(p.xz)-args.x,p.y);
    return length(q)-args.y;
}

//Cylinder: uncapped cylinder at pos c.xy and radius c.z
float fCylinder( vec3 p, vec3 c ){
    return length(p.xz-c.xy)-c.z;
}

//Cone: a cone
float fCone( vec3 p, vec3 args ){
    // c must be normalized
    float q = length(p.xy);
    return dot(args.xy,vec2(q,p.z));
}

//Plane: a plane of normal n
float fPlane( vec3 p, vec3 n ){
    // n must be normalized
    return dot(p,n.xyz);
}

//HexPrism: a hexagonal prism of radius args.x and width args.y
float fHexPrism( vec3 p, vec3 args ){
    vec3 q = abs(p);
    return max(q.z-args.y,max((q.x*0.866025+q.y*0.5),q.y)-args.x);
}

//TriPrism: a triangular prism of size args.x and width args.y
float sdTriPrism( vec3 p, vec3 args ){
    vec3 q = abs(p);
    return max(q.z-args.y,max(q.x*0.866025+p.y*0.5,-p.y)-args.x*0.5);
}

// Capsule: A Cylinder with round caps on both sides
float fCapsule(vec3 p, float r, float c) {
    return mix(length(p.xz) - r, length(vec3(p.x, abs(p.y) - c, p.z)) - r, step(c, abs(p.y)));
}

//CappedCylinder: a cylinder of radius args.x and length args.y
float fCappedCylinder( vec3 p, vec3 args )
{
    vec2 d = abs(vec2(length(p.xz),p.y)) - args.xy;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

//CappedCone: a capped cone of
float fCappedCone( vec3 p, vec3 c )
{
    vec2 q = vec2( length(p.xz), p.y );
    vec2 v = vec2( c.z*c.y/c.x, -c.z );
    vec2 w = v - q;
    vec2 vv = vec2( dot(v,v), v.x*v.x );
    vec2 qv = vec2( dot(v,w), v.x*w.x );
    vec2 d = max(qv,0.0)*qv/vv;
    return sqrt( dot(w,w) - max(d.x,d.y) ) * sign(max(q.y*v.x-q.x*v.y,w.y));
}

//Ellipsoid: and ellipsoid
float fEllipsoid( vec3 p, vec3 r ){
    return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z);
}

////////////////////////////////////////////////////////////////
//
//                DOMAIN MANIPULATION OPERATORS
//
////////////////////////////////////////////////////////////////
//
// Conventions:
//
// Everything that modifies the domain is named pSomething.
//
// Many operate only on a subset of the three dimensions. For those,
// you must choose the dimensions that you want manipulated
// by supplying e.g. <p.x> or <p.zx>
//
// <inout p> is always the first argument and modified in place.
//
// Many of the operators partition space into cells. An identifier
// or cell index is returned, if possible. This return value is
// intended to be optionally used e.g. as a random seed to change
// parameters of the distance functions inside the cells.
//
// Unless stated otherwise, for cell index 0, <p> is unchanged and cells
// are centered on the origin so objects don't have to be moved to fit.
//
//
////////////////////////////////////////////////////////////////



// Rotate around a coordinate axis (i.e. in a plane perpendicular to that axis) by angle <args.x>.
// Read like this: R(p.xz, a) rotates "x towards z".
// This is fast if <args.x> is a compile-time constant and slower (but still practical) if not.
void pR(inout vec2 p, vec3 args) {
    p = cos(args.x)*p + sin(args.x)*vec2(p.y, -p.x);
}

// Shortcut for 45-degrees rotation
void pR45(inout vec2 p) {
    p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

// Repeat space along one axis. Use like this to repeat along the x axis:
// <float cell = pMod1(p.x,5);> - using the return value is optional.
float pMod1(inout float p, vec3 args) {
    float size = args.x;
    float halfsize = size*0.5;
    float c = floor((p + halfsize)/size);
    p = mod(p + halfsize, size) - halfsize;
    return c;
}


// Repeat around the origin by a fixed angle.
// For easier use, num of repetitions is use to specify the angle.
float pModPolar(inout vec2 p, vec3 args) {
    float repetitions = args.x;
    float angle = 2.*PI/repetitions;
    float a = atan(p.y, p.x) + angle/2.;
    float r = length(p);
    float c = floor(a/angle);
    a = mod(a,angle) - angle/2.;
    p = vec2(cos(a), sin(a))*r;
    // For an odd number of repetitions, fix cell index of the cell in -x direction
    // (cell index would be e.g. -5 and 5 in the two halves of the cell):
    if (abs(c) >= (repetitions/2.)) c = abs(c);
    return c;
}

// Repeat in two dimensions
vec2 pMod2(inout vec2 p, vec3 args) {
    vec2 size = args.xy;
    vec2 c = floor((p + size*0.5)/size);
    p = mod(p + size*0.5,size) - size*0.5;
    return c;
}

// Repeat in three dimensions
vec3 pMod3(inout vec3 p, vec3 size) {
    vec3 c = floor((p + size*0.5)/size);
    p = mod(p + size*0.5, size) - size*0.5;
    return c;
}


////////////////////////////////////////////////////////////////
//
//             OBJECT COMBINATION OPERATORS
//
////////////////////////////////////////////////////////////////
//
// We usually need the following boolean operators to combine two objects:
// Union: OR(a,b)
// Intersection: AND(a,b)
// Difference: AND(a,!b)
// (a and b being the distances to the objects).
//
// The trivial implementations are min(a,b) for union, max(a,b) for intersection
// and max(a,-b) for difference. To combine objects in more interesting ways to
// produce rounded edges, chamfers, stairs, etc. instead of plain sharp edges we
// can use combination operators. It is common to use some kind of "smooth minimum"
// instead of min(), but we don't like that because it does not preserve Lipschitz
// continuity in many cases.
//
// Naming convention: since they return a distance, they are called fOpSomething.
// The different flavours usually implement all the boolean operators above
// and are called fOpUnionRound, fOpIntersectionRound, etc.
//
// The basic idea: Assume the object surfaces intersect at a right angle. The two
// distances <a> and <b> constitute a new local two-dimensional coordinate system
// with the actual intersection as the origin. In this coordinate system, we can
// evaluate any 2D distance function we want in order to shape the edge.
//
// The operators below are just those that we found useful or interesting and should
// be seen as examples. There are infinitely more possible operators.
//
// They are designed to actually produce correct distances or distance bounds, unlike
// popular "smooth minimum" operators, on the condition that the gradients of the two
// SDFs are at right angles. When they are off by more than 30 degrees or so, the
// Lipschitz condition will no longer hold (i.e. you might get artifacts). The worst
// case is parallel surfaces that are close to each other.
//
// Most have a float argument <r> to specify the radius of the feature they represent.
// This should be much smaller than the object size.
//
// Some of them have checks like "if ((-a < r) && (-b < r))" that restrict
// their influence (and computation cost) to a certain area. You might
// want to lift that restriction or enforce it. We have left it as comments
// in some cases.
//
// usage example:
//
// float fTwoBoxes(vec3 p) {
//   float box0 = fBox(p, vec3(1));
//   float box1 = fBox(p-vec3(1), vec3(1));
//   return fOpUnionChamfer(box0, box1, 0.2);
// }
//
////////////////////////////////////////////////////////////////

//Union: joins two primitives
vec2 fOpUnion(vec2 a, vec2 b){
    if (a.x < b.x){
        return a;
    }
    return b;
}

//Intersect: takes only where both primitives intersect
vec2 fOpIntersect(vec2 a, vec2 b){
    if (a.x > b.x){
        return a;
    }
    return b;
}

//Difference: subtract a from b;
vec2 fOpDifference(vec2 a, vec2 b){
    return vec2(fOpIntersect(a, vec2(-b.x, b.y)).x,a.y);
}

// The "Chamfer" flavour makes a 45-degree chamfered edge (the diagonal of a square of size <args.x>):
vec2 fOpUnionChamfer(vec2 a, vec2 b, vec3 args) {
    float r = args.x;
    vec2 u = fOpUnion(a,b);
    float m = u.x;
    //if ((a < r) && (b < r)) {
    return vec2(min(m, (a.x - r + b.x)*sqrt(0.5)),u.y);
    //} else {
    return u;
    //}
}

// Intersection has to deal with what is normally the inside of the resulting object
// when using union, which we normally don't care about too much. Thus, intersection
// implementations sometimes differ from union implementations.
vec2 fOpIntersectionChamfer(vec2 a, vec2 b, vec3 args) {
    float r = args.x;
    vec2 i = fOpIntersect(a,b);
    float m = i.x;
    if (r <= 0.) return i;
    if (((-a.x < r) && (-b.x < r)) || (m < 0.)) {
        return vec2(max(m, (a.x + r + b.x)*sqrt(0.5)),i.y);
    } else {
        return i;
    }
}

// Difference can be built from Intersection or Union:
vec2 fOpDifferenceChamfer (vec2 a, vec2 b, vec3 args) {
    return fOpIntersectionChamfer(a, vec2(-b.x,b.y), args);
}

// The "Round" variant uses a quarter-circle to join the two objects smoothly:
vec2 fOpUnionRound(vec2 a, vec2 b, vec3 args) {
    float r = args.x;
    vec2 u = fOpUnion(a,b);
    float m = u.x;
    if ((a.x < r) && (b.x < r) ) {
        return vec2(min(m, r - sqrt((r-a.x)*(r-a.x) + (r-b.x)*(r-b.x))),u.y);
    } else {
        return u;
    }
}

vec2 fOpIntersectionRound(vec2 a, vec2 b, vec3 args) {
    float r = args.x;
    vec2 i = fOpIntersect(a,b);
    float m = i.x;
    if ((-a.x < r) && (-b.x < r)) {
        return vec2(max(m, -(r - sqrt((r+a.x)*(r+a.x) + (r+b.x)*(r+b.x)))),i.y);
    } else {
        return i;
    }
}

vec2 fOpDifferenceRound (vec2 a, vec2 b, vec3 args) {
    return fOpIntersectionRound(a, vec2(-b.x,b.y), args);
}


////////////////////////////////////////////////////////////////
//
//                         MISC
//
////////////////////////////////////////////////////////////////

float noise(vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.-2.*f);
    
    float n = p.x + p.y*157. + 113.*p.z;
    
    vec4 v1 = fract(753.5453123*sin(n + vec4(0., 1., 157., 158.)));
    vec4 v2 = fract(753.5453123*sin(n + vec4(113., 114., 270., 271.)));
    vec4 v3 = mix(v1, v2, f.z);
    vec2 v4 = mix(v3.xy, v3.zw, f.y);
    return mix(v4.x, v4.y, f.x);
}

vec2 rand(vec2 p)
{
    vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);
    
}

vec3 tex3DWithFilter(sampler2D tex, vec3 pos, in vec3 nor)
{
    vec3 uvw = pos/10*8;
    
    // calc texture sampling footprint
    vec3 ddx_uvw = uvw + dFdx( uvw );
    vec3 ddy_uvw = uvw + dFdy( uvw );
    
    int MaxSamples = 10;
    int sx = 1 + int( clamp( 4.0*length(ddx_uvw-uvw), 0.0, float(MaxSamples-1) ) );
    int sy = 1 + int( clamp( 4.0*length(ddy_uvw-uvw), 0.0, float(MaxSamples-1) ) );
    
    vec3 no = vec3(0.0);
    
    for( int j=0; j<MaxSamples; j++ )
        for( int i=0; i<MaxSamples; i++ )
        {
            if( j<sy && i<sx )
            {
                vec2 st = vec2( float(i), float(j) ) / vec2( float(sx),float(sy) );
                no += tex3D( tex, uvw + st.x*(ddx_uvw-uvw) + st.y*(ddy_uvw-uvw), nor );
            }
        }
    
    return no / float(sx*sy);
}

vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n ){
    
    n = max((abs(n) - 0.2)*7., 0.001);
    n /= (n.x + n.y + n.z );
    return (texture(tex, p.yz)*n.x + texture(tex, p.zx)*n.y + texture(tex, p.xy)*n.z).xyz;
}

vec3 texBump( sampler2D tx, in vec3 p, in vec3 n, float bf){
    
    const vec2 e = vec2(0.002, 0);
    
    mat3 m = mat3( tex3DWithFilter(tx, p - e.xyy, n), tex3D(tx, p - e.yxy, n), tex3D(tx, p - e.yyx, n));
    
    vec3 g = vec3(0.299, 0.587, 0.114)*m;
    g = (g - dot(tex3DWithFilter(tx,  p , n), vec3(0.299, 0.587, 0.114)) )/e.x; g -= n*dot(n, g);
    
    return normalize( n + g*bf );
    
}


////////////////////////////////////////////////////////////////
//
//                 OBJECT MODIFICATION OPERATORS
//
////////////////////////////////////////////////////////////////



vec2 fBumpify(vec3 p, vec2 shape, vec3 args){
    float bumpiness = args.x;
    if (shape.x > 1) return shape;
    mat3 M = mat3(0.28862355854826727, 0.6997227302779844, 0.6535170557707412,
                  0.06997493955670424, 0.6653237235314099, -0.7432683571499161,
                  -0.9548821651308448, 0.26025457467376617, 0.14306504491456504);
    vec3 p1 = M*p;
    vec3 p2 = M*p1;
    float n1 = noise(p1*5.);
    float n2 = noise(p2*10.);
    float n3 = noise(p1*20.);
    float n4 = noise(p1*40.);
    float rocky = 0.1*n1*n1 + 0.05*n2*n2 + 0.02*n3*n3 + 0.01*n4*n4;
    return vec2(shape.x + rocky*bumpiness, shape.y);
}

////////////////////////////////////////////////////////////////
// The end of HG_SDF library
////////////////////////////////////////////////////////////////
