#include "Shapes.h"


void Plane::test(Ray& ray, HitData& hit)
{
	float t = (this->d - this->n.Dot(ray.o)) / this->n.Dot(ray.d);
	if ((hit.t > t || hit.t < 0) && t >= 0)
	{
		hit.t = t;
		hit.color = this->c;
		hit.lastNormal = this->n;
		hit.lastShape = this;
	}

}
Vec Plane::normal(Vec &point)
{
	return this->n;
}
Plane::Plane(Vec normal, float _d, Color color)
{
	this->n = normal;
	// Maybe not necessary
	this->n.Normalize();
	this->d = _d;
	this->c = color;
}
//-------------------------------------------------------

void Sphere::test(Ray& ray, HitData& hit)
{
	float b = ray.d.Dot(ray.o - this->center);
	float c = (ray.o - this->center).Dot(ray.o - this->center) - this->radius2;
	if (b*b >= c)
	{
		float t = -b - std::sqrtf(b*b-c);
		if (hit.t > t || hit.t < 0)
		{
			hit.t = t;
			hit.color = this->c;
			hit.lastNormal = this->normal(ray.o+ray.d*t);
			hit.lastShape = this;
		}
	}

}
Vec Sphere::normal(Vec &point)
{
	Vec n = point - this->center;
	n.Normalize();
	return n;
}
Sphere::Sphere(Vec _center, float _radius, Color _color)
{
	this->center = _center;
	this->radius = _radius;
	this->radius2 = radius*radius;
	this->c = _color;
}

//------------------------------------------------------

float Triangle::det(const Vec & u, const Vec & v, const Vec & w) const
{
	return
			u.x*(v.y*w.z - v.z*w.y)
		-	u.y*(v.x*w.z - v.z*w.x)
		+	u.z*(v.x*w.y - v.y*w.x);
}

void Triangle::test(Ray & ray, HitData & hit)
{
	Vec s = ray.o - this->p1;
	Vec nd = Vec(0.f,0.f,0.f) - ray.d;
	float scl = 1 / this->det(nd, this->edge0, this->edge1);
	float u = scl*this->det(nd, s, this->edge1);
	if (u >= 0)
	{
		float v = scl*this->det(nd, this->edge0, s);
		if (u + v <= 1.f && v >= 0)
		{
			float t = scl*this->det(s, this->edge0, this->edge1);
			if (hit.t > t || hit.t < 0)
			{
				hit.t = t;
				hit.color = this->c;
				hit.lastNormal = this->nor;
				hit.lastShape = this;
			}
		}
	}
	
}

Triangle::Triangle(Vec _p1, Vec _p2, Vec _p3, Color _color)
{
	this->p1 = _p1;
	this->p2 = _p2;
	this->p3 = _p3;
	this->edge0 = _p2 - _p1;
	this->edge1 = _p3 - _p1;
	this->nor = Vec(
		this->edge0.y*this->edge1.z - this->edge0.z*this->edge1.y,
		this->edge0.z*this->edge1.x - this->edge0.x*this->edge1.z,
		this->edge0.x*this->edge1.y - this->edge0.y*this->edge1.x);
	this->nor.Normalize();
	this->c = _color;
}

void OBB::test(Ray & ray, HitData & hit)
{
	float tmin = -INFINITY;
	float tmax = INFINITY;
	Vec p = this->Bcenter - ray.o;
	Vec *base[3] = { &this->Bu,&this->Bv,&this->Bw };
	float half[3] = { this->halfBu,this->halfBv,this->halfBw };
	bool reject = false;
	for (int i = 0; i < 3 && !reject; i++)
	{
		float e = base[i]->Dot(p);
		float f = base[i]->Dot(ray.d);
		if (std::fabs(f) > 0.f)
		{
			float t1 = (e + half[i]) / f;
			float t2 = (e - half[i]) / f;
			if (t1 > t2)
			{
				float temp = t1;
				t1 = t2;
				t2 = temp;
			}
			tmin = std::fmaxf(tmin, t1);
			tmax = std::fminf(tmax, t2);
			if (tmin > tmax || tmax < 0)
				reject = true;
		}
		else if (-e - half[i] > 0 || -e + half[i] < 0)
			reject = true;
	}
	if (!reject)
	{
		float t;
		if (tmin > 0)
			t = tmin;
		else
			t = tmax;
		//if (hit.lastShape == nullptr || t < hit.t)
		if (hit.t > t || hit.t < 0)
		{
			hit.t = t;
			hit.color = this->c;
			hit.lastNormal = this->normal(ray.o + ray.d*hit.t);
			hit.lastShape = this;
		}
	}
}

Vec OBB::normal(Vec & point)
{
	Vec nor;
	if (std::fabs(this->Bu.Dot(point) - this->Bu.Dot(this->Pu)) < 0.001f)
		nor = this->Bu;
	else if (std::fabs(this->Bu.Dot(point) - this->Bu.Dot(this->Puo)) < 0.001f)
		nor = Vec(0, 0, 0) - this->Bu;
	else if (std::fabs(this->Bv.Dot(point) - this->Bv.Dot(this->Pv)) < 0.001f)
		nor = this->Bv;
	else if (std::fabs(this->Bv.Dot(point) - this->Bv.Dot(this->Pvo)) < 0.001f)
		nor = Vec(0, 0, 0) - this->Bv;
	else if (std::fabs(this->Bw.Dot(point) - this->Bw.Dot(this->Pw)) < 0.001f)
		nor = this->Bw;
	else if (std::fabs(this->Bw.Dot(point) - this->Bw.Dot(this->Pwo)) < 0.001f)
		nor = Vec(0, 0, 0) - this->Bw;

	return nor;
}

OBB::OBB(Vec b, Vec b1, Vec b2, Vec b3, float Hu, float Hv, float Hw, Color _color)
{
	this->Bcenter = b;
	this->Bu = b1;
	this->Bv = b2;
	this->Bw = b3;
	this->Bu.Normalize();
	this->Bv.Normalize();
	this->Bw.Normalize();
	this->halfBu = Hu;
	this->halfBv = Hv;
	this->halfBw = Hw;
	this->Pu =	b + this->Bu*Hu;
	this->Puo = b - this->Bu*Hu;
	this->Pv =	b + this->Bv*Hv;
	this->Pvo = b - this->Bv*Hv;
	this->Pw =	b + this->Bw*Hw;
	this->Pwo = b - this->Bw*Hw;

	this->c = _color;
}

OBB::OBB(Vec b, float Hu, float Hv, float Hw, Color _color)
{
	this->Bcenter = b;
	this->halfBu = Hu;
	this->halfBv = Hv;
	this->halfBw = Hw;
	this->Bu = Vec(1, 0, 0);
	this->Bv = Vec(0, 1, 0);
	this->Bw = Vec(0, 0, 1);
	this->c = _color;
}

Color Shape::shade(Vec & light, const Vec & cam, Ray & r, HitData & h)
{
	Vec ambient = Vec(50, 50, 50);
	Vec point = (r.o + r.d*h.t);
	point = light - point;
	point.Normalize();
	float angle = h.lastNormal.Dot(point);
	if (angle < 0.f)
		angle = 0.f;

	BYTE red =		std::fminf(h.color.r*angle + h.color.r*(ambient.z/255.f), 255);
	BYTE green =	std::fminf(h.color.g*angle + h.color.g*(ambient.x/255.f), 255);
	BYTE blue =		std::fminf(h.color.b*angle + h.color.b*(ambient.y/255.f), 255);
	
	return Color(
		(BYTE)red,
		(BYTE)green,
		(BYTE)blue);
}
