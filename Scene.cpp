//--------------------------------------------------------------------------------------
// File: Lab01-ImplementHere.cpp
//
// DV1541/DV1542 Lab01: RayTracing
//
// Copyright (c) Stefan Petersson 2010. All rights reserved.
//--------------------------------------------------------------------------------------
#include "stdafx.h"
#include "ToScreen.h"
#include <vector>
#include <algorithm>
#include "Shapes.h"
#include <assert.h> 

#pragma comment(lib, "shapes.lib")

extern ToScreen* g_ToScreen;
#define TORAD(X) ((X)*PI/180.0f)

// container for each shape in the scene
std::vector<Shape*> shapes;



// adds all the shapes to a vector, that is used in the Update function.
void buildScene()
{
	if (shapes.size() == 0)
	{
		// TESTING
		/*Triangle test = Triangle({ 2, 0, 0 }, { 4,0,0 }, { 3, 0, 2 }, { 200,0,0 });
		Ray testray = Ray(Vec(3, 10, 1), Vec(0, -1, 0));
		test.test(testray, HitData());*/

		// a plane facing the camera, passing through point (0,0,1000)
		shapes.push_back(dynamic_cast<Shape*>(new Plane(Vec(0, 0, -1), -1000, Color(0, 30, 30))));

		// a sphere at (400,400,300), with radius 200
		shapes.push_back(dynamic_cast<Shape*>(new Sphere(Vec(400, 400, 300), 200, Color(100, 100, 0))));

		// two triangles 
		shapes.push_back(dynamic_cast<Shape*>(new Triangle({ 350, 100, 200 }, { 300,100,200 }, { 400,700, 30 }, { 200,0,0 })));
		shapes.push_back(dynamic_cast<Shape*>(new Triangle({ 100,300,0 }, { 150,300,0 }, { 100,100,0 }, { 0,0,255 })));
		shapes.push_back(dynamic_cast<Shape*>(new Triangle({ 150,300,0 }, { 150,100,0 }, { 100,100,0 }, { 0,255,255 })));


		float size = 20;
		Vec p0{ -1, 0, 0 }; Vec p1{ 1, 0, 0 }; Vec p2{ 0,2*size,0 };
		for (unsigned int i = 1; i <= 20; i++)
		{
			Vec off { size*2 * i, 700.0, 0.0 };
			p0.x = cosf(TORAD(i * (90/20.0f))) * size; p0.z = -sinf(TORAD(i * (90/20.0f))) * size;
			p1.x = -p0.x; p1.z = -p0.z;
			shapes.push_back(dynamic_cast<Shape*>(new Triangle(p0 + off, p1 + off, p2 + off, { 255,255,0 })));
		}
		
		// one OBB touching the sphere on the side
		// base
		Vec b1{ 1, 0, 0 };
		Vec b2{ 0, 1, 0 };
		Vec b3{ 0, 0, 1 };
		float angle{0.5};
		// rotate around Z the basis
		b1.x = cosf(angle); b1.y = -sinf(angle);
		b2.x = sinf(angle); b2.y = cosf(angle);
		shapes.push_back(dynamic_cast<Shape*>(new OBB(Vec(g_ToScreen->mScreenWidth/2, g_ToScreen->mScreenHeight/2, 100), b1, b2, b3, 50, 50, 50, { 0,255,0 })));
		// further rotate around X the basis
		float tempY = b1.y * cosf(angle) + b1.z * sinf(angle);
		b1.z = b1.y * -sinf(angle) + b1.z * cosf(angle);
		b1.y = tempY;
		tempY = b2.y * cosf(angle) + b2.z * sinf(angle);
		b2.z = b2.y * -sinf(angle) + b2.z * cosf(angle);
		b2.y = tempY;
		shapes.push_back(dynamic_cast<Shape*>(new OBB(Vec(200,600,400), b1,b2,b3, 100, 100, 100, {255,0,0})));
	}
}


// This method should not be modified for the purpose of the assignment.
void Update(float deltaTime)
{
	buildScene();
	// randomise the scene.
	// we do this to verify that no matter in what order the primitives are tested,
	// they should always return the same intersection points
	std::random_shuffle(shapes.begin(), shapes.end()); 

	// THE CODE BELOW CANNOT BE CHANGED FOR THIS ASSIGNMENT
	//loop through every screen pixel from bottom to top
	for(unsigned int y = 0; y < g_ToScreen->mScreenHeight; y++)
	{
		for(unsigned int x = 0; x < g_ToScreen->mScreenWidth; x++)
		{
			// Black color by default
			Color c = { 0, 0, 0 };

			// hitData will be reset every pixel.
			// hitData is sent to intersection tests to return intersection info
			// and ALWAYS carries the closest intersection for each pixel
			HitData hitData;

			// Create ray for current pixel, in ORTHOGRAPHIC fasion
			// Ray from (x,y,-10) with direction (0,0,1)
			Ray r(Vec((float)x, (float)y, -10), Vec(0,0,1));
			
			for (auto& shapeIt : shapes)
				shapeIt->test(r, hitData);

			// if we hit anything, we shade it
			if (hitData.lastShape != nullptr)
			{
				Shape* s = dynamic_cast<Shape*>(hitData.lastShape);
				// light source position at the center of the screen, z=-1000
				// camera position at the center of the screen, z = -10, 
				// normal at surface is queried from the surface and hitData
				c = s->shade(Vec(g_ToScreen->mScreenWidth/2, g_ToScreen->mScreenHeight/2, -1000), Vec(g_ToScreen->mScreenWidth/2, g_ToScreen->mScreenHeight/2, -10), r, hitData);
			}
			// hitData should contain the closest intersection point, 
			// and of course the colour of the shape
			// If there was no intersection, the colour will be BLACK
			g_ToScreen->SetPixelColor( 
				x, 
				(g_ToScreen->mScreenHeight-1-y), 
				c.r, c.g, c.b );
		}
	}
	//update texture surface with new pixel data
	g_ToScreen->Update(deltaTime);
}