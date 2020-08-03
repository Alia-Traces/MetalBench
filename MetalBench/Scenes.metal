	//
	//  Scenes.metal
	//  MetalBench
	//
	//  Created by Alia on 31/07/2020.
	//

#include <metal_stdlib>
using namespace metal;

#define GLASSSCENE		0
#define FORESTSCENE	1

constant int SCENEFUNC [[function_constant(0)]];
//constant int SCENEFUNC = FORESTSCENE;
constant int INCREMENTALRENDER = SCENEFUNC == FORESTSCENE;
constant int REALTIMERENDER = !INCREMENTALRENDER;

#define GroundMaterial	{0, {.3,.8,.3}, .2, .2, 0, 0, .8}

	// Camera values
	//#define CameraPosition float3(fract(ray.origin.w/8.) * -30 + 15,12,-25)
	//#define CameraTarget 	float3(fract(ray.origin.w/8.) * -40 + 20,4,0)
	//#define CameraPosition float3(-11,40,1)
	//#define CameraPosition float3(sin(ray.origin.w)*15, 12, cos(ray.origin.w)*15)

	//#define Zoom 4//2 //(uv.x<0?128:2) // Lens zoom.
	//#define Fisheye 1.7 // Use with wider lenses (zoom of 1 or less)
	//#define Aperture 0 //(0.5*(sin(ray.origin.w)*.5+.5)) // For depth of field
	//#define SoftLens 2
	//#define BokehShape .25 // Changes the bokeh shape between gaussian (1), disk (0.25) and ring (0)

	// The scene can have a main light (a sphere)
	//#define LightPosition float3(-500,200,-400)*5
	//#define LightSize 500 // Sphere radius
	//#define LightCol float3(150,120,100)*4

	// Materials
#define GlassClear	 	0
#define GlassBlue	 	1
#define GlassPink		2
#define Leaf			3
#define MatTest		4

#define GroundMat	{0,{.27,.16,.1},.1} //	{0, {1, .8, .6}, .2}

	// 2d inplace rotation. p = float2 to rotate, a = angle
#define Rotate2D(p,a) p=cos(a)*p+sin(a)*float2(-p.y,p.x);
#define mod(a, b) (a - b * floor(a/b))

	// Minor values
#undef Epsilon
#define Epsilon 1e-3 // epsilon

#pragma mark - Structs

struct Ray { // Ray
	float4 origin;
	float3 dir;
	bool inside; // Used to indicate ray is inside surface (e.g. glass)
};

struct Box {
	float3 a, b;
};

struct Triangle {
	float3 v, a, b;
};

struct Material{
	int id; // 0: standard, 1: light, 2: glass
	float3 color;
		// Standard materials have 2 layers, allowing for e.g. matte base and reflective coat
	float baseSmoothness, coatSmoothness, coatOpacity, fresnel, coatSaturation;
};

	/// Stores the result of an intersection test - normal, distance, material
struct Intersection {
	float4 result; // Result: xyz = normal, w = dist
	Material material; // Material
			   //	S scene; // S number
};


#pragma mark - Hashing

	/// Basic hash function
float3 hash(float3 p){
	p=fract(p*float3(443.897,441.423,437.195));
	p+=dot(p,p.yxz+19.19);
	return fract((p.xxy+p.yxx)*p.zyx);
}

#pragma mark - Sampling

	/// Returns a point on a sphere, r is in 0..1 range
float3 pointOnSphere(float2 r) {
	r=float2(6.283185*r.x,2*r.y-1);
	return float3(sqrt(1+Epsilon-r.y*r.y)*float2(cos(r.x),sin(r.x)),r.y); // 1.001 required to avoid NaN
}

	/// Returns a cosine weighted sample
float3 lambertSample(float3 n,float2 r) {
	return normalize(n*(1+Epsilon)+pointOnSphere(r)); // 1.001 required to avoid NaN
}

#pragma mark Prototypes

Intersection intersectScene(thread Ray &ray, float3 k);

#pragma mark Misc

#pragma mark - Lighting

bool lightSample(thread Ray &ray, Intersection hit, float3 k, float4 i, thread float3 &o) {
	float3 lightPosition, lightSize;
	if (SCENEFUNC == FORESTSCENE) {
		lightPosition = float3(-280,250,200)*2;
		lightSize = 50;
	}
	
		// Get r to random point on light
	float3 d = normalize(lightPosition + pointOnSphere(k.xy) * k.z * lightSize - ray.origin.xyz);
	
		// Get cosine of r and normal
	float c=dot(d,hit.result.xyz);
	
		// Reject if facing away from light
	if(c<=0) return 0;
	
		// See if the ray hits a light (can be any light or the sky)
	Ray ray2 = ray;
	ray2.dir = d;
	hit = intersectScene(ray2,k);
	
	if(hit.material.id==1){
		o+=i.xyz*c*hit.material.color*(hit.material.color/(4*M_PI_F*hit.result.w));
	}
	if(hit.material.id==2){
		ray=ray2;
		return 1;
	}
	return 0;
}

#pragma mark - Materials

	/// Applies glass, including fresnel reflection
void applyGlass(thread Ray&ray, Intersection hit, float3 k, thread float4 &rayColouration){
	
		// Fresnel term
	float fresnel = 1 - max(0., -dot(ray.dir, hit.result.xyz));
	fresnel = pow(fresnel, 2.);
	
		// Randomly reflect or refract, probability based on fresnel term
		// If non-fresnel, refract
	if (k.z > fresnel) {
			// refraction
			// Index of refraction
		float ior = ray.inside ? 1.5 : 1./1.5;
		
			// Find refraction angle, accounting for surface roughness
		float3 rayDir = normalize( // Ray dir must be normalised
					  mix( // Mix between...
					      refract(ray.dir, hit.result.xyz, ior), // the refracted angle
					      lambertSample(-hit.result.xyz, k.xy), // and a random angle projected into the surface
					      hit.material.baseSmoothness // based on how smooth the glass surface is
					      )
					  );
		
			// Test for total internal reflection
		if(dot(hit.result.xyz, rayDir) < 0) {
				// Not TIR, we're OK to refract
				// Step through suface along normal
			ray.origin.xyz -= hit.result.xyz * Epsilon * 4;
			
				// Set the ray direction
			ray.dir=rayDir;
			
				// Flip the inside value as we pass through the surface
			ray.inside=!ray.inside;
			return;
		}
	}
	
		// Ray failed to refract, therefore reflection
		// Step away from the surface along the normal
	ray.origin.xyz+=hit.result.xyz*Epsilon*2;
	
		// Standard reflection
	ray.dir=reflect(ray.dir,hit.result.xyz);//mix(reflect(ray.dir,hit.result.xyz),lambertSample(hit.result.xyz,k.xy),hit.material.baseSmoothness);
}


	/// This is hacky, I need to rework it, but it works for basic materials for now
void applyMaterial(thread Ray&ray, thread Intersection &hit, float3 k, thread float4 &rayColouration, thread float3 &outputColour){
	if (SCENEFUNC == GLASSSCENE) {
			// useROI, useMainLight, angular (0=volume)
			// b, cs, o, f, s
			// There's no refraction so step away from suRayFactorace to prevent re-intersection
		ray.origin.xyz+=hit.result.xyz*Epsilon*2;
		
			// Fresnel value
		float f=1+dot(ray.dir,hit.result.xyz),c; // 0 head on, 1 at oblique angle
		f*=f*f*hit.material.fresnel; // Fresnel is dependent on the coating
		
			// Add the coat opacity. This means the outer suRayFactorace will be visible, but still respects fresnel
		f+=(hit.material.coatOpacity);
		
			// Split randomly between coat and base
			//	Should write this without conditionals - coatLayer = 0 | 1, then just change base values to base or coat accordingly
		c=step(k.z,f);
		
			// Colour by m, blended into whe (no colouration) according to coat colouration (allows for coloured reflections)
		rayColouration.xyz*=mix(hit.material.color,mix(1,hit.material.color,hit.material.coatSaturation),c);
		
			// We'll use base m only, so set that to base or coat
		hit.material.baseSmoothness=pow(mix(hit.material.baseSmoothness,hit.material.coatSmoothness,c),.5);
			//	hit.material.baseSmoothness = 0.0;
			// If there's a main light, sample it
			//	if(lightSample(ray, hit, k, rayColouration, outputColour)) return;
			// ROI ray?
			//	if (hit.material.baseSmoothness < k.z) {
			////		if (roiCast(ray, hit.result.xyz, k)) return;
			//	}
			//	lightSample(r,h,k,i,o);
		
			//		ray.dir = LightPosition + pointOnSphere(k.xy) * k.z * LightSize;
			//		ray.dir = normalize(ray.dir - ray.origin.xyz);
			//		return;
			//	}
		
			//
			// Standard suRayFactorace handling
			// Properties: base roughness, coat roughness, coat opacity, coat fresnel
		
			// Get random direction
		float3 s=pointOnSphere(k.xy),n=normalize(hit.result.xyz+s*.9*(1-pow(hit.material.baseSmoothness,.125)));
		ray.dir=normalize(reflect(ray.dir,n)*hit.material.baseSmoothness // mirror
				  +s*pow(k.z+Epsilon,hit.material.baseSmoothness*hit.material.baseSmoothness*10) // rough
				  );
		
			// If normal points away from r, flip it
		if(dot(ray.dir,hit.result.xyz)<0)ray.dir=normalize(hit.result.xyz*(1+Epsilon)+s*pow(k.z+Epsilon,0.25));
	} else {
			// There's no refraction so step away from suRayFactorace to prevent re-intersection
		ray.origin.xyz+=hit.result.xyz*Epsilon*2;
		
		
			// Fresnel value
		float f=1+dot(ray.dir,hit.result.xyz); // 0 head on, 1 at oblique angle
		f*=f*f*hit.material.fresnel; // Fresnel is dependent on the coating
		
			//	float spec = k.z < hit.material.coatOpacity ? hit.material.coatSmoothness : hit.material.baseSmoothness;
		float spec = mix(hit.material.baseSmoothness, hit.material.coatSmoothness,  step(1-hit.material.coatOpacity, k.z));
		spec = f > k.z ? hit.material.coatSmoothness : spec;
		spec = spec * 20 + .5;
		
			//	rayColouration.xyz *= hit.material.color;
		rayColouration.xyz *= pow(hit.material.color, 1/spec);
		
		if(lightSample(ray, hit, k, rayColouration, outputColour)) return;
			//	return;
		k = pointOnSphere(k.xy);
		float3 n = normalize(k + hit.result.xyz*spec);
			//	float d = dot(n, hit.result.xyz);
		
			// If normal points in same direction as ray, colourise and flip
		if (dot(n, hit.result.xyz)<0) {
			n = -n;
		}
			//	d = dot(n, ray.dir);
		if (dot(n, ray.dir)>0){
			rayColouration.xyz *= pow(hit.material.color, 1/spec);
			ray.dir = normalize(hit.result.xyz + k);
		} else ray.dir = reflect(ray.dir, n);
	}
}

void check(Ray ray, thread Intersection &hit, float s) {
	ray.origin.xyz = ray.origin.xyz + ray.dir.xyz * hit.result.w;
	ray.origin.xz = mod(floor(ray.origin.xz * s), 2.0);
	hit.material.color *= floor(fmod(ray.origin.x + ray.origin.z, 2.0) * .95 + .05);
}
	//
	//void polka(Ray ray, thread Intersection &hit, float s){ // C+O
	//	ray.origin.xyz = ray.origin.xyz + ray.dir.xyz * hit.result.w;
	//	hit.material.color *= step(0.35, length(mod(ray.origin.xz * s, 1.0) - 0.5));
	//}
	//
	//void lightPolka(Ray ray, thread Intersection &hit, float s, float ls){ // C+O
	//	ray.origin.xyz = ray.origin.xyz + ray.dir.xyz * hit.result.w;
	//	hit.material.color *= (1.-step(ls, length(mod(ray.origin.xz * s, 1.0) - 0.5)))*2.;
	////	float2 p = abs(fract(ray.origin.xz * s)-.5);
	////	p = step(ls, p);
	////	hit.material.color = float3(1-min(p.x, p.y));
	//	hit.material.id = 1;
	//}

#pragma mark - returnay intersection functions

	/// Ground plane intersection
void intersectGround(Ray ray, thread Intersection &hit, Material material){
	ray.origin.w = -ray.origin.y / ray.dir.y;
	if(ray.origin.w > 0 && ray.origin.w < hit.result.w){
		hit.result = float4(0,1,0, ray.origin.w);
		hit.material = material;
	}
}

	/// Sphere intersection
void intersectSphere(Ray ray, float4 sphere, thread Intersection &hit, Material material){
	ray.origin.xyz -= sphere.xyz;
	ray.origin.w = dot(ray.dir.xyz, ray.origin.xyz) * 2;
	float a = dot(ray.origin.xyz, ray.origin.xyz) - sphere.w * sphere.w;
	a = ray.origin.w * ray.origin.w - 4 * a;
	if (a < 0) return;
	a = sqrt(a);
	float2 g = (float2(-a, a) - ray.origin.w) / 2;
	a = g.x < 0 ? g.y : g.x;
	sphere.w *= sign(g.x);
	if (a> hit.result.w || a < 0) return;
	hit.result = float4((ray.dir.xyz * a + ray.origin.xyz) / sphere.w, a);
	hit.material = material;
}

	/// Sphere intersection points
	//float2 intersectSphere(float3 o, float3 d, float4 sphere, thread float3 &n){
	//	o -= sphere.xyz;
	//	float w = dot(d, o) * 2,
	//	a = dot(o, o) - sphere.w * sphere.w;
	//	a = w * w - 4 * a;
	//	if (a < 0) return -MAXFLOAT;
	//	a = sqrt(a);
	//	float2 g = (float2(-a, a) - w) / 2;
	//	n = normalize(d * (g.x < 0 ? g.y : g.x) + o);// * sign(g.x);
	//	return g;
	//}

	/// Cube intersection
	//void intersectCube(Ray r,Box c,thread Intersection&hit,Material m){
	//	float3 a=(c.a-r.origin.xyz)/r.dir, // near
	//	b=(c.b-r.origin.xyz)/r.dir, // far
	//	f=max(a,b), // furthest
	//	n=min(a,b); // nearest
	//	float x=min(f.x,min(f.y,f.z)), // furthest plane
	//	d=max(n.x,max(n.y,n.z)), // nearest plane
	//	o=d<0?x:d; // nearest in front
	//	if(isnan(n.x)|d>=x|o>hit.result.w|o<0)return; // d>=x = invalid, o>t = behind other geometry, o<0 behind
	//	hit.result.w=o;
	//	hit.result.xyz=normalize(step(Epsilon,abs(a-hit.result.w))-step(Epsilon,abs(b-hit.result.w)))*sign(d);
	//	hit.material=m;
	//}


	/// Intersection test for triangle
	//void intersectTriangle(Ray r, Triangle t, thread Intersection &hit, Material m){
	//	float3 p=cross(r.dir,t.b),q,s;
	//	float e=dot(t.a,p),u,v;
	//	if(e<Epsilon)return;
	//
	//	float f=1/e;
	//
	//	s=r.origin.xyz-t.v;
	//	u=dot(s,p)*f;
	//	float i=step(0,u)*(1-step(1,u));
	//
	//	q=cross(s,t.a);
	//	v=dot(r.dir,q)*f;
	//	i*=step(0,v)*(1-step(1,u+v));
	//	if(i==0)return;
	//
	//	u=dot(t.b,q)*f;
	//
	//	f=step(0,-u);
	//	u=u*(1-f)+(f*MAXFLOAT);
	//	if(u>hit.result.w)return;
	//	p=normalize(cross(t.a,t.b));
	//	hit = {float4(p*sign(e),u),m};
	//}

void intersectPlane(Ray r, float4 plane, thread Intersection &hit, Material m) {
	float dist = dot(plane.xyz * plane.w - r.origin.xyz, plane.xyz) / dot(r.dir, plane.xyz);
	if (dist < hit.result.w && dist > 0) {
		hit.result = float4(plane.xyz, dist);
		hit.material = m;
	}
}

void sphereSliceMirror(Ray ray, thread Intersection &hit, Material mat, float4 sphere, float4 plane, float planeOffset) {
	float3 p = ray.origin.xyz - sphere.xyz;
		// Intersect sphere
	float w = dot(ray.dir, p) * 2, w2,
	a = dot(p, p) - sphere.w * sphere.w;
	a = w * w - 4 * a;
	if (a < 0) return;
	a = sqrt(a);
	float2 g = (float2(-a, a) - w) / 2;
	if (max(g.x,g.y)<0) return;
	a = g.x<0?g.y:g.x; // Distance to first intersection
	if (a<0) return;
	float3 n = normalize(ray.dir * a + p) * sign(g.x);
	
	w = dot(plane.xyz * (plane.w + planeOffset) - p, plane.xyz) / dot(ray.dir, plane.xyz); // Dist to plane
	w2 = dot(-plane.xyz * (plane.w - planeOffset) - p, -plane.xyz) / dot(ray.dir, -plane.xyz); // Dist to plane
	
	bool facingPlane = dot(ray.dir, plane.xyz) < 0;//, facingPlane2 = !facingPlane;
	float3 pn = plane.xyz * (facingPlane * 2 - 1), pn2 = -pn;
	
	if (g.x>0) {
		if (max(w, w2) < 0) {
				// past far plane
			return;
		} else if (a > max(w, w2)) {
				// Ray intersects sphere past far plane
			return;
		} else {
			a = max(a,min(w, w2));
			n = a==w ? pn : (a==w2 ? pn2 : n);
		}
	} else {
		if (max(w, w2) < 0) {
				// Past furthest point
			return;
		} else if (min(w, w2) < 0) {
				// Inside slice, take minimum to sphere, further plane
			a = min(a, max(w, w2));
			n = a==w ? pn : (a==w2 ? pn2 : n);
		} else if (a>w) {
			return;
		}
	}
	
	if (length(p + ray.dir * a) > sphere.w + Epsilon) return; // Exit if ray intersects sphere past the cutoff
	if (a>hit.result.w) return; // Exit if beyond nearest surface
	if (dot(ray.dir, n) > 0) {n = -n;}
	
		// Set material
	hit.result = float4(n, a);
	hit.material = mat;
}

#pragma mark - Distance functions

	/// Box, origin is o, size is s
float boxDist(float3 p, float3 o, float3 s){
	float3 q = abs(p-o) - s/2;
	return length(max(q,0.)) + min(max(q.x,max(q.y,q.z)),0.);
}

	/// Torus
float torus(float3 p, float2 t){
	float2 q = float2(length(p.xz)-t.x,p.y);
	return length(q)-t.y;
}

float cylinderDist(float3 p, float3 o, float2 s){
	p -= o;
	float2 d = abs(float2(length(p.xz),p.y))-s;
	return min(max(d.x,d.y),0.) + length(max(d,0.));
}

float smin(float a, float b, float k) {
	float res = exp2(-k*a) + exp2(-k*b);
	return -log2(res)/k;
}

#pragma mark Forest scene

	// Dist, stencil
float2 leaf(float3 p, float3 o, float3 s) {
		// Twist and bend
	p -= o;
	Rotate2D(p.xz, s.z * 6);
	p.z -= s.z * 5;
	Rotate2D(p.yz, .5*p.z/20 * s.y);
	p.z += s.z * 5;
	Rotate2D(p.xy, sin(p.z/3) * s.y*.3);
	float3 q = p;
	
	float e = s.x*3+3,
	d = max(length(p.xy) - (p.z+e)/50, length(p)- 10+s.y*3),f; // d = stalk
								   //	d = max(length(p.xy) - sqrt((max(0.,p.z+e)))/10, length(p)- 10+s.y*3),f; // d = stalk
	
	p.x = abs(p.x); // Mirror
	p.z = fract(p.z+p.x)-.5; // Ribs
	f = length(p.yz) - .15; // f = rib tubes
	
	p = q;
	p.x = abs(p.x) + e; // Oval leaf shape
	
		//	e = length(p) - 10+s.y*3; // e = leaf volume (sphere-ish)
	e = max(length(p.xz) - 10+s.y*3, -p.y); // e = leaf volume (Tube-ish)
	p.y = abs(p.y)-.1; // Leaf thickness
	f = max(f, e); // Clamp rib tubes to leaf outline
		       //	e = max(e, p.y); //Leaf is now flat
	d = min(d, min(max(e, p.y), f)); // Add leaf to ribs and stalk
	return {d, e};
}

	// Dist, stencil, colour
float3 leafLayer(float3 p) {
	p.xz /= 20;
	float3 s = hash(floor(p.xzz)+.1),
	l = float3(leaf(float3(fract(p.xz) *20 - 10, p.y).xzy, {0,.1,0}, s), -(s.x+s.y+s.z+1));
	p-=.5;
	Rotate2D(p.xz, 3.15);
	s = hash(floor(p.xzz)+.1);
	s = float3(leaf(float3(fract(p.xz) *20 - 10, p.y).xzy, {0,.1,0}, s), -(s.x+s.y+s.z+1));
		//	p.xz = fract(p.xz) *20 - 10;
	return l.x<s.x?l:s;
}

float2 ring(float3 p, float3 o) {
	p -= o;
	return {torus(p*float3(4,1.6,4), {12,1})/4, 1};
}

#pragma mark - scene

	/// Distance function
float2 df(float3 p, float3 dir) {
	float3 q=p, e;
	float2 d = MAXFLOAT;
		//	return matTest(p);
	p.xz += 100;
	p.y -= .6;
		//	p.xz = mod(p.xz, 20) - 10;
		//	float2 d = leaf(p, {0,.1,0}, {5, 10}),e;
	float h = 0;
	for (int i=0; i<5; i++){ //5
		e = leafLayer(p);
			// Stencil anything ABOVE the leaf!
		d.x = max(d.x, -e.y);
		d = d.x<e.x?d:e.xz;
		p += {7,-.5,7};
		Rotate2D(p.xz, 1.);
	}
	e.xy = ring(q, float3(12,2,56) + float3(-4, 1.1, -5));
	d = d.x<e.x?d:e.xy;
	return d;
}

void rayMarch(Ray ray, thread Intersection &hit, float3 k){
		// Current position
	float3 p = ray.origin.xyz ;//+ ray.dir * k.z*0 ;
	
		// How much to step the ray by. Inverted when ray is inside
	float scale = 1;// ray.inside?-1:1;
			//	bool inside = boundingDist(p);
			//	if (!inside) return;
			//	scale *= .5;
			//	scale*=1; // If you need shorter steps (inaccurate distance function) set it here
			//
			//
			// Marching loop
	for(int i=0; i<200; i++) {
		if (p.y > 7 && ray.dir.y > 0) return;
			// Check if we're still in bounds. If not, either trace to next bounds or exit
			//		bool insideNow = boundingDist(p, ray.dir) < .1;
			//		insideNow = true;
			//		if (!insideNow) {
			//			// Trace to bounds
			//			ray.origin.xyz = p;
			//			Intersection h = {float4(MAXFLOAT), {}};
			//			float b = intersectBounds(ray, h);
			//			if (b >= hit.result.w) break; // Ray didn't intersect bounds before next surface
			//			p += ray.dir * (b + Epsilon); // trace to next bounds
			//			continue;
			//		}
			// Total distance ray travelled (used for early termination and returning intersection distance)
		float totalDist = length(ray.origin.xyz-p);// + k.z*0;
		
			// Early termination if ray travelled further than the last intersection, or if out of bounds
		if(totalDist>hit.result.w) return;
			//		if(p.y>TileSize * 18 & ray.dir.y>0) return;
			//		if(oob(p, ray.dir)) return;
			//		if(totalDist>hit.result.w || abs(p.y) > dfThickness+Epsilon) return;
		
			// Get the distance (x = dist, y = material ID)
		float2 dist=df(p, ray.dir);
		
			// Check if we hit a surface
		if (abs(dist.x) < Epsilon * 2) {
				// intersection
				// Get normals
			float2 e = float2(Epsilon * .1, 0.);
			float3 n=normalize(
					   float3(
						  df(p+e.xyy, ray.dir).x-df(p-e.xyy, ray.dir).x,
						  df(p+e.yxy, ray.dir).x-df(p-e.yxy, ray.dir).x,
						  df(p+e.yyx, ray.dir).x-df(p-e.yyx, ray.dir).x
						  )
					   );
				// Restore if nans appear
			if(any(isnan(n))) continue;
			
				//			p.xz = fract(p.xz / 1);
			float3 c = floor(p.xzz/10)+10;//, q = hash(c);
			if (dist.y<0) {
					//				c = mix(float3(1,.5,.6), float3(.6, .7, .4), (abs(dist.y)-1) / 3);
				c = hash(dist.yyy);
				c = mix(float3(1,.5,.2), float3(.8, .9, .2), c.x);
				dist.y = Leaf;
			}
			
				//			c = {1,.6,.2};
			
				// An array of materials
			Material m[]={
				GroundMat, // Ground
				{0,{1,.8,0},1,0,0,0}, // Ring
						      //				{0, 1, 1},
				{2,{0,0,0},0}, // Diamond
				{0, c, .02}, // Leaf
					     //				{0,hash(c)*.5+.5, fract(c.x / 4.), fract(x), q.z, 1, 1}, //4 MatTest
				{0,{1,.6,.2}}, //4 MatTest
					       //				{0,1,0}, //5 Pride flag
					       //				GroundMaterial, // ground
					       //				{0,{1,.1,.1},0.2,.8,.2,0,0.5}, // RedBrick
					       //				{0,{.1,1,.1},0.2,.8,.2,0,0.5}, //7 GreenBrick
					       //				{0, {1,1,1}, .5}, // 8 MarbleMat
					       //				{2,{0},1}, // 9 sea
					       //				{0,1,.3}, // 10 White
					       //				{1,{2,0.5,.5}}, // 11 pink light
					       //				{2, {1,0,0},0.2} // 12 red glass
			};
			
				// Set the material
			hit.material = m[int(dist.y)];
			
				// Flip normals if inside object
			if (hit.material.id == 2 && ray.inside) n = -n;
			
				// Set the intersection result
			hit.result=float4(n,totalDist);
			
				//			hit.material = {1, i / 30};
			return;
		}
		
			// Step along ray
		p+=ray.dir*dist.x*scale;
	}
}

Intersection intersectScene(thread Ray &ray, float3 k){
	Intersection hit, testHit;
	float3 n, o, p, blue, pink;
	float thickness;
	bool l;
	
	switch (SCENEFUNC) {
		case GLASSSCENE:
				// Create an intersection result. This is the background, so no intersection needed. The normal is
				// the opposite of the ray dir and distance is 'far away'. Material is undefined, we'll set that later.
			hit = {float4(-ray.dir, MAXFLOAT), {1,mix(float3(1,.5,.3)*2, float3(0,.0,1), pow(max(0., ray.dir.y), 0.4))}};
			
			intersectGround(ray, hit, {0,1}); // White
			
				// Back wall
			intersectPlane(ray, {0,0,-1,-10}, hit, {0,1,.3});
			
				// Intersect marble
			n = normalize(float3(-1,3,-1)), o={17,5.01,0},
			blue = {.08,.05,0},
			pink = {0,.06,.04};
			thickness = 1-Epsilon * 8;
			
				// Jiggled marble
			n = normalize(float3(1,0,0));
			sphereSliceMirror(ray, hit, {2,blue},	float4(o+n*4,5), 	float4(n,thickness), -4);
			sphereSliceMirror(ray, hit, {2,pink},	float4(o+n*6,5), 	float4(n,thickness), -2);
			sphereSliceMirror(ray, hit, {2,0}, 		float4(o+n*2,5), 	float4(n,thickness), 0);
			sphereSliceMirror(ray, hit, {2,pink},	float4(o-n*4,5), 	float4(n,thickness), 2);
			sphereSliceMirror(ray, hit, {2,blue},	float4(o-n*8,5), 	float4(n,thickness), 4);
			
				//	 Sliced marble
			n = normalize(float3(-1,-.0,-.55));
			sphereSliceMirror(ray, hit, {2,blue},	float4(4,5.01,0,5), 		float4(n,1), -4);
			sphereSliceMirror(ray, hit, {2,pink},	float4(2,5.01,0,5), 		float4(n,1), -2);
			sphereSliceMirror(ray, hit, {2}, 		float4(0,5.01,0,5), 		float4(n,1), 0);
			sphereSliceMirror(ray, hit, {2,pink},	float4(-2,5.01,0,5), 	float4(n,1), 2);
			sphereSliceMirror(ray, hit, {2,blue},	float4(-4,5.01,0,5), 	float4(n,1), 4);
			
				// Joined marble
			n = normalize(float3(3,1,.5));
			sphereSliceMirror(ray, hit, {2,blue},	float4(-17,5.01,0,5), 	float4(n,thickness), -4);
			sphereSliceMirror(ray, hit, {2,pink},	float4(-17,5.01,0,5), 	float4(n,thickness), -2);
			sphereSliceMirror(ray, hit, {2,0}, 		float4(-17,5.01,0,5), 	float4(n,thickness), 0);
			sphereSliceMirror(ray, hit, {2,pink},	float4(-17,5.01,0,5), 	float4(n,thickness), 2);
			sphereSliceMirror(ray, hit, {2,blue},	float4(-17,5.01,0,5), 	float4(n,thickness), 4);
			
				// Lights
			
				//	intersectSphere(ray, float4(-40, 40, 0, 10), hit, {1, 15});// blue
			intersectSphere(ray, float4(45, 20, 6, 4), hit, {1, {40,35,30}});// blue
			intersectSphere(ray, float4(-20, 30, 5, 5), hit, {1, {10,15,20}});// blue
			break;
			
		case FORESTSCENE:
			
			hit = {float4(-ray.dir, MAXFLOAT), {1,mix(float3(.8,.8,.5), float3(1,.7,.3)*2, pow(max(0., ray.dir.y), 0.5))/4}};
			
				// Ground
			intersectGround(ray, hit, {0,{.27,.16,.1},.1});
			
				// Intersect sun mask
			float3 lightPosition = float3(-280,250,200)*2;
			
			testHit = hit;
			intersectPlane(ray, float4(normalize(lightPosition), length(lightPosition) / 3), testHit, {});
			if (testHit.result.w < hit.result.w) {
				p = fract((ray.origin.xyz + ray.dir * testHit.result.w+150)/30)-.5;
				l = length(p)-.45 < 0;
				p = fract((ray.origin.xyz + ray.dir * testHit.result.w * 1.2-10)/5)-.5;
				l &= length(p)-.5 < 0;
				if (l) {
						// Intersect sun
					intersectSphere(ray, float4(lightPosition, 50), hit, {1, float3(150,70,20)*4});
				}
			}
			
				// Raymarch (do this last for better performance)
			rayMarch(ray, hit, k);
			break;
	}
	
	return hit;
}

#pragma mark - Tracing
	//================//

	//float focalLength(Ray ray, float3 camPos, float3 camTarget) {
	//	float3 n =normalize(camPos - camTarget);
	//	return -dot(ray.origin.xyz - camPos, n) / dot(ray.dir, n);
	//}

	/// Sets the camera up. uv = screen position, k = random value
Ray setupCamera(Ray ray, float2 uv, float3 k){
	float3 camPos, camTarget;
	
	float zoom = 2, aperture = 0, softLens = 0, fisheye = 0, bokehShape = 0.25;
	
	float t = ray.origin.w / 12., segment = floor(mod(t, 3.));
	t = fract(t);
	
	switch (SCENEFUNC) {
		case GLASSSCENE:
				//			segment = 2;
			if (segment == 0) {
					// Pan
				camPos = float3(t * -30 + 15,12,-25);
				camTarget = float3(t * -40 + 20,4,0);
				zoom = 4;
				softLens = 2;
				fisheye = 1.7;
			} else if (segment == 1) {
					// Focus pull
				camPos = float3(-40, 25, 4);
				camTarget = float3(t * -35 + 20,4,0);
				zoom = 6 - fract(t) * 4;
				aperture = 4.0;
			} else {
					// Rotation
				camPos = float3(0, 7, -7);
				camTarget = float3(-18, 5-t*5, 0);
				Rotate2D(camPos.xz, -t*2 + 1);
				camPos.x -= 15;
				zoom = 1;
				fisheye = 1.0;
			}
			break;
		case FORESTSCENE:
			camPos = float3(-20, 33, 30);
			camTarget = float3(12,2,56);
			zoom = 2.5;
			fisheye = 1.5;
			aperture = 1.0;
			softLens = 1.5;
			break;
	}
	
		// Basic lens zoom first
	uv /= zoom;
	
		// get the ray dir
	ray.dir = camTarget - camPos;
	
		// f = focal length
	float fl = length(ray.dir), a;
	ray.dir = normalize(ray.dir);
	
		// This transforms k into a random point in a sphere
		//	k = hash(k);
	k = pointOnSphere(k.xy) * pow(k.z, bokehShape) * (aperture + length(uv*zoom) * softLens);
	
		// Add random sphere point * aperture to camera position for DoF
		//	ray.origin.xyz = CameraPosition + k * Aperture;// * length(uv); // Can uncomment this to create soft focus at edges only
	
		// Update the ray direction, then project back from the camera target to the camera plane
		//	ray.dir=normalize(CameraTarget-ray.origin.xyz);
		//	ray.origin.xyz=CameraTarget-ray.dir*f;
	ray.origin.xyz = camPos;
		//	ray.origin.xyz += k * pow(1-abs(dot(normalize(k), ray.dir)), 2);
	ray.dir = normalize(camTarget - ray.origin.xyz);
		// Transform the camera to account for uv
	a = rsqrt(1 - ray.dir.y * ray.dir.y);
	float3x3 c = {
		float3(-ray.dir.z, 0, ray.dir.x) * a,
		float3(-ray.dir.x * ray.dir.y, 1 - ray.dir.y * ray.dir.y, -ray.dir.y * ray.dir.z) * a,
		-ray.dir
	};
	a=length(uv);
	
		// Scaling for fisheye distortion
	float f = a * fisheye;
	
		// Last bit of uv transform
	ray.dir=normalize(c*mix(float3(uv,-1),float3(uv/a*sin(f),-cos(f)),.4));
	
		//	fl = focalLength(ray);
	float3 t2 = normalize(camPos - camTarget); // plane facing camera
	fl = -dot(ray.origin.xyz - camTarget, t2) / dot(ray.dir, t2);
	t2 = ray.origin.xyz + ray.dir * fl;
	ray.origin.xyz += k;
	ray.dir = normalize(t2 - ray.origin.xyz);
	
	return ray;
}

float3 traceRay(float2 uv, float pixelSize, float t, int rayCount) {
	Ray ray;
	
	int bounces;
	switch (SCENEFUNC) {
		case GLASSSCENE:
			bounces = 15;
			break;
		case FORESTSCENE:
			bounces = 4;
			break;
	}
	
	float3 lightSum=0, // Light sum, we add lights to this and return it at the end
	glassCol,
	k=hash(float3(uv + t,t)); // initial random value
				  //	k=hash(float3(t *0+ 0.1)); // initial random value
				  //	k = hash(float3(tgp, t+.1));
				  //	k = hash(uv.yyy);
	
		// Go through rays
	for(int j=0; j<rayCount; j++){
		ray.inside = 0; // Set to 1 if the camera is inside a glass object!
		
			// This gives us 2 values of stratified sampling plus one random value on z
			//		k.xy = float2(
			//			      fmod(float(j), rayFactor),
			//			      floor(float(j) / rayFactor)
			//			      ) / rayFactor + k.xy / rayFactor;
		
			// To save memory and registers, use append the time to the ray origin
		ray.origin.w=t + (float(j)/float(rayCount*20));
		
			// setup camera, get ray
			// We add a small random value to uv here, this gives us anti-aliasing
		ray = setupCamera(ray, uv + hash(k).xy * pixelSize * 2, k);
		
			// This stores the ray colouration. If the ray hits a surface, this is multiplied by the surface colour.
			// That stores the light absorbed by each surface the ray intersects
			// If we hit a light, we just multiply the light colour by the ray colouration to get a light value for the whole path
			// The w value stores hue if we do spectral rendering
		float4 rayColouration = {1,1,1};
		
			// iterate through Bounces
		for(int bounce=0; bounce<bounces; bounce++) {
				// intersect and move ray
			Intersection hit = intersectScene(ray, k);
			
				// If the ray is inside an object, subtract the material colour scaled by the distance the ray just travelled
				// This gives us accurate coloured glass, including when the ray bounces multiple times inside the object
			if(ray.inside) {
				rayColouration.xyz = max(0., rayColouration.xyz - (glassCol * hit.result.w));
			}
			
				// Move along ray to surface
			ray.origin.xyz += ray.dir * hit.result.w;
			
				// Apply the material
				// Light
			if(hit.material.id == 1) {
					// Add light and terminate the ray
					//				float3 c = hit.material.color.xyz, l=LightCol;
					//				if (c.x+c.y+c.z == l.x+l.y+l.z) { break; }
				lightSum += hit.material.color * rayColouration.xyz;
				if (bounce == 0) lightSum += lightSum * ((k.z * 2 - 1) * (2./(rayCount)));
				break;
			}
			
				// Glass
			if(hit.material.id == 2) {
				glassCol = hit.material.color;
				applyGlass(ray, hit, k, rayColouration);
			}
			
				// Standard material
			if(hit.material.id == 0) {
				applyMaterial(ray, hit, k, rayColouration, lightSum);
			}
			
				// Early exit if h light or r nearly expired
				// This improves performance in scenes with dark surfaces, as the ray gets terminated if it's too dark
				//This gets rid of speckles!!!! WHY
				//			if(all(rayColouration.xyz < .3)) { break;}
				//			if(any(rayColouration.xyz < 0.1)) {break;}
				//			if (any(isnan(ray.dir))) { lightSum += {100,0,100}; break; }
				//			if (any(lightSum < 0)) { lightSum += {1,0,1}; }
				//
				// Sometimes useful for debug...
				//			if (Bounces >= 2.0) break;
				//			k = hash(floor(ray.origin.xyz*1000));
		}
		
			// New hash value for next ray
		k=hash(k);
	}
	return lightSum;
}

#pragma mark YUV
float3 rgbToYUV709 (float3 rgb) {
	return float3x3(
			float3(.2126, .7152, .00722),
			float3(-.1145, -.3854, .5),
			float3(.5, -.4541, -.04584)
			) * rgb;
}

float3 yuv709ToRGB (float3 yuv) {
	return float3x3(
			float3(1.069, .1289, 1.5749),
			float3(1.069, -0.058, -0.468),
			float3(1.069, 1.984, 0)
			) * yuv;
}

#pragma mark - Compute functions
	//=========================//
	//
kernel void k(
	      texture2d<float, access::write> o [[texture(0)]],
	      texture2d<float, access::read> i [[texture(1)]],
	      constant float&t[[buffer(0)]],
	      constant int &rays[[buffer(1)]],
	      constant uint2 &gridOffset[[buffer(2), function_constant(INCREMENTALRENDER)]],
	      uint2 g[[thread_position_in_grid]]
	      ){
	
	float brightness = 1,
	gamma = 1./2.2,
	contrast = 1.0,
	contrastBias = 0.5,
	blackPoint = 0.0,
	whitePoint = 1.0,
	warmth = 1,
	saturation = 1.0;
	
	float2 shadowTint = {0,0},
	highlightTint = {0,0},
	tint = {0,0};
	
	switch (SCENEFUNC) {
		case GLASSSCENE:
			gamma = .6;
			contrastBias = 0.4;
			saturation = 1.2;
			shadowTint = {-.1,.4};
			highlightTint = {-.1,-.8};
			tint = {-.2,-.05};
			break;
		case FORESTSCENE:
			whitePoint = 2.0;
			gamma = .5;
			contrast = 1.2;
			contrastBias = 0.3;
			saturation = 0.5;
			warmth = 0.7;
			tint = {-0.3,0.1};
			shadowTint = {-0.3,-0.2};
			highlightTint = {1,1};
			break;
	}
	float2 r=float2(o.get_width(),o.get_height());
	float2 u;
	if (INCREMENTALRENDER) {
		u = (float2((g + gridOffset)*2)-r)/float2(r.y, -r.y);
	} else {
		u = (float2((g)*2)-r)/float2(r.y, -r.y);
	}
	
	float3 p = traceRay(
			    u,
			    5 / r.x, // AA
			    fmod(t, M_PI_F * 20.0),
			    rays
			    );
	p = p*pow(saturate(2.2-length(u)),.75) / rays;
	p = pow(p, gamma);
	
	p = rgbToYUV709(p);
	
	p.z *= brightness;
	p.z = (p.z-contrastBias) * contrast + contrastBias;
	p.z = saturate(p.z / whitePoint) / (1 - blackPoint) + blackPoint;
	
	p.xy = pow(abs(p.xy),warmth) * sign(p.xy) * saturation;
	
	p.xy += mix(shadowTint, highlightTint, pow(max(0., p.z), 1.3)) + tint;
	
	p = yuv709ToRGB(p);
	
	if (SCENEFUNC == FORESTSCENE) {
		float4 l = i.read(g + gridOffset);
		l += float4(p, 1);
		o.write(l, g + gridOffset);
	} else {
		// Standard output
		o.write(float4(p, 1), g);
	}
	
}

kernel void displayIntegrated (
			       texture2d<float,access::read> i [[texture(0)]],
			       texture2d<float,access::write> o [[texture(1)]],
			       uint2 g[[thread_position_in_grid]]
			       ){
	float4 p = i.read(g);
	p.xyz /= p.a;
	p.a = 1;
	o.write(p, g);
}
