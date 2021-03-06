#include "Brdfread.hpp"
#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{

    // Change the definition here to change resolution
    Scene scene(300, 300);

    Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Material* light = new Material(DIFFUSE, 
    (80.0f * Vector3f(0.747f+0.058f, 0.747f+0.258f, 0.747f) +
     150.6f * Vector3f(0.740f+0.287f,0.740f+0.160f,0.740f) +
      180.4f *Vector3f(0.737f+0.642f,0.737f+0.159f,0.737f)));
    light->Kd = Vector3f(5.65f);


    Material* measured = new Material(MEASURED, Vector3f(0.0f));
    //use Kd as measured brdf
    const char *filename = argv[1];
	  double* brdf;
    if (!read_brdf(filename, brdf)) 
	  {
		  fprintf(stderr, "Error reading %s\n", filename);
		  exit(1);
	  }
    measured->brdf_p=brdf;
	// read brdf
	  
    measured->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    //MeshTriangle floor("../models/cornellbox/floor.obj", white);
    //MeshTriangle shortbox("../models/cornellbox/shortbox.obj", white);
    //MeshTriangle tallbox("../models/cornellbox/tallbox.obj", white);
    //MeshTriangle left("../models/cornellbox/left.obj", red);
    //MeshTriangle right("../models/cornellbox/right.obj", green);
    MeshTriangle light_("../models/cornellbox/triangle.obj", light);
    MeshTriangle bunny("../models/bunny/bunny1.obj", measured);

    //scene.Add(&floor);
    //scene.Add(&shortbox);
    //scene.Add(&tallbox);
    //scene.Add(&left);
    //scene.Add(&right);
    scene.Add(&light_);
    scene.Add(&bunny);

    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}