//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    auto p=this->intersect(ray);
    auto wo=-ray.direction;
  
    if(p.happened==false)return Vector3f(0);

    if(p.m->hasEmission())return p.m->getEmission();
    Intersection Inter;
    float pdf;
    sampleLight(Inter,pdf);

    auto xx=Inter.coords;
    auto NN=Inter.normal;
    auto newRay=Ray(p.coords,normalize(xx-p.coords));
    auto ws=newRay.direction;

    auto testRayInter=intersect(newRay);
    Vector3f L_dir,L_indir;
    if(fabs(testRayInter.distance-(xx-p.coords).norm())<=EPSILON*10)
    {
        L_dir=Inter.emit *p.m->eval(wo,ws,p.normal)*dotProduct(ws,p.normal)*dotProduct(-ws,NN)/powf((xx-p.coords).norm(),2)/pdf;
    }

    float num=get_random_float();
    if(num<RussianRoulette)
    {
        auto wi=p.m->sample(wo,p.normal);
        auto rr=Ray(p.coords,normalize(wi));
        auto inte=intersect(rr);
        if(inte.happened==false)return Vector3f(0);
        if(!inte.m->hasEmission())
        {
            L_indir=castRay(rr,depth+1)*p.m->eval(wo,wi,p.normal)*dotProduct(wi,p.normal)/p.m->pdf(wo,wi,p.normal)/RussianRoulette;
        }
    }
    return L_indir+L_dir;//; 
}