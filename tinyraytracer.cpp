#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"
#include "material.h"
#include "model.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


int envmap_width, envmap_height;
std::vector<Vec3f> envmap;

const int nbFrame = 360;

struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float &r, const Material &m) : center(c), radius(r), material(m) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

struct Box {
    Vec3f point1;
    Vec3f point2;
    Material material;

    Box(const Vec3f &p1, const Vec3f &p2, const Material &m) : point1(p1), point2(p2), material(m) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        /*
        Vec3f tMin = Vec3f((point1.x - orig.x) / dir.x, (point1.y - orig.y) / dir.y, (point1.z - orig.z) / dir.z);
        Vec3f tMax = Vec3f((point2.x - orig.x) / dir.x, (point2.y - orig.y) / dir.y, (point2.z - orig.z) / dir.z);
        Vec3f t1;
        Vec3f t2;
        if((tMin.x + tMin.y + tMin.z < tMax.x + tMax.y + tMax.z)){
            t1 = tMin;
            t2 = tMax;
        }
        else{
            t1 = tMax;
            t2 = tMin;
        }

        t0 = std::max(std::max(t1.x, t1.y), t1.z);
        float tFar = std::min(std::min(t2.x, t2.y), t2.z);
        if(t0>tFar)return false;
        return true;

        */

        //Vec3f t1 = std::min(tMin, tMax);
        //Vec3f t2 = std::max(tMin, tMax);
        //t0 = std::max(std::max(t1.x, t1.y), t1.z);
        //float tFar = std::min(std::min(t2.x, t2.y), t2.z);
        //if(t0>tFar)return false;
        //return true;

        float tmin = (point1.x - orig.x) / dir.x; 
        float tmax = (point2.x - orig.x) / dir.x; 
    
        if (tmin > tmax) std::swap(tmin, tmax); 
    
        float tymin = (point1.y - orig.y) / dir.y; 
        float tymax = (point2.y - orig.y) / dir.y; 
    
        if (tymin > tymax) std::swap(tymin, tymax); 
    
        if ((tmin > tymax) || (tymin > tmax)) 
            return false; 
    
        if (tymin > tmin) 
            tmin = tymin; 
    
        if (tymax < tmax) 
            tmax = tymax; 
    
        float tzmin = (point1.z - orig.z) / dir.z;  
        float tzmax = (point2.z - orig.z) / dir.z; 
    
        if (tzmin > tzmax) std::swap(tzmin, tzmax); 
    
        if ((tmin > tzmax) || (tzmin > tmax)) 
            return false; 
    
        if (tzmin > tmin) 
            tmin = tzmin; 
    
        if (tzmax < tmax) 
            tmax = tzmax; 

    
        t0 = tmin;
        if (t0 < 0) return false;
        return true; 
        
    }
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

Vec3f refract(const Vec3f& I, const Vec3f& N, const float eta_t, const float eta_i = 1.f) { // Snell's law
    float cosi = -std::max(-1.f, std::min(1.f, I * N));
    if (cosi < 0) return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(1, 0, 0) : I * eta + N * (eta * cosi - sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, Model *model, const std::vector<Box> &boxes, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }
    

    float boxes_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < boxes.size(); i++) {
        float dist_i;
        if (boxes[i].ray_intersect(orig, dir, dist_i) && dist_i < boxes_dist && dist_i < spheres_dist) {
            boxes_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - boxes[i].point2).normalize(); // a verif
            material = boxes[i].material;
        }
    }

    //model
    float models_dist = std::numeric_limits<float>::max();

    if(model->intersect_bbox(orig, dir)){
        for (int i=0; i<model->nfaces(); i++) {
            float dist_i;
            if (model->ray_triangle_intersect(i,orig, dir, dist_i) && dist_i < models_dist && dist_i < boxes_dist && dist_i < spheres_dist) {
                models_dist = dist_i;
                hit = orig + dir*dist_i;
                Vec3f point1 = model->point(model->vert(i, 0));
                Vec3f point2 = model->point(model->vert(i, 1));
                Vec3f point3 = model->point(model->vert(i, 2));
                N = (cross(point2-point1, point3-point1)).normalize();
                material = model->getMaterial();
            }
        }
    }

    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(dir.y)>1e-3)  {
        float d = -(orig.y+4)/dir.y; // the checkerboard plane has equation y = -4
        Vec3f pt = orig + dir*d;
        if (d>0 && fabs(pt.x)<10 && pt.z<-10 && pt.z>-30 && d < models_dist && d < boxes_dist && d < spheres_dist) {
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0,1,0);
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(.3, .3, .3) : Vec3f(.3, .2, .1);
        }
    }
    return std::min(spheres_dist, std::min(boxes_dist, std::min(models_dist, checkerboard_dist))) < 1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, Model *model, const std::vector<Box> &boxes, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth=0) {
    Vec3f point, N;
    Material material;

    if (depth>4 || !scene_intersect(orig, dir, model, boxes, spheres, point, N, material)) {
        //return Vec3f(0.2, 0.7, 0.8); // background color
/*
        int i = ((((dir.x/(width/(float)height)/tan(fov/2.))/2)*(float)width  - 1))-0.5;
        int j = -((dir.y/tan(fov/2.)/2)*(float)height - 1)-0.5;
        
        std::cout << "i: "<< i << std::endl;
        std::cout << "j: " << j << std::endl;

        */
        
        
        int i = ((int)((atan2(dir.z, dir.x) / (2 * M_PI) + 0.5) * envmap_width));
        int j = (int)(acos(dir.y) / M_PI * envmap_height);

        return envmap[(i+j*envmap_width)%(envmap_width*envmap_height)];
        
/*
        Vec3f p = dir;
        float theta = acosf(p.y / p.norm());
        float phi = atan2f(p.z, p.x) + M_PI;

        int i = theta / (M_PI) * (envmap_height);
        int j = phi / (2 * M_PI) * (envmap_width);

        //return envmap[std::min(i + j * envmap_width, envmap_width*envmap_height-1)];
        return envmap[(i + j * envmap_width)%(envmap_width*envmap_height)];

*/

        /*
        Vec3f p = dir;
        float theta = acosf(p.y / p.norm());
        float phi = atan2f(p.z, p.x) + M_PI;

        int y = theta / (M_PI) * (envmap_height);
        int x = phi / (2 * M_PI) * (envmap_width);

        return envmap[x + y * envmap_width];
        */

    }

    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // offset the original point to avoid occlusion by the object itself
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, model, boxes, spheres, lights, depth + 1);
    Vec3f refract_color = cast_ray(refract_orig, refract_dir, model, boxes, spheres, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, model, boxes, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflect_color*material.albedo[2] + refract_color*material.albedo[3];
}

float frame = 0;

float depthMat = 2000.f;

Vec3f m2v(Matrix m) {
    return Vec3f(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}

Matrix v2m(Vec3f v) {
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depthMat/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depthMat/2.f;
    return m;
}

Matrix lookat(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye-center).normalize();
    Vec3f x = (up^z).normalize();
    Vec3f y = (z^x).normalize();
    Matrix res = Matrix::identity(4);
    for (int i=0; i<3; i++) {
        res[0][i] = x[i];
        res[1][i] = y[i];
        res[2][i] = z[i];
        res[i][3] = -center[i];
    }
    return res;
}

Matrix rot_z(const float angle) {
    Matrix R = Matrix::identity(4);
    R[0][0] = R[2][2] = cos(angle);
    R[2][0] = sin(angle);
    R[0][2] = -R[2][0];
    return R;
}

Matrix rot_y(const float angle) {
    Matrix R = Matrix::identity(4);
    R[0][0] = R[1][1] = cos(angle);
    R[1][0] = sin(angle);
    R[0][1] = -R[1][0];
    return R;
}

Matrix rot_x(const float angle) {
    Matrix R = Matrix::identity(4);
    R[1][1] = R[2][2] = cos(angle);
    R[1][2] = sin(angle);
    R[2][1] = -R[1][2];
    return R;
}

void render(Model *model, const std::vector<Box> &boxes, const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    const int width    = 1024;
    const int height   = 768;
    const int fov      = M_PI/2.;
    std::vector<Vec3f> framebuffer(width*height);

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            //float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
            //float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            float x =  (i + 0.5) -  width/2.;
            float y = -(j + 0.5) + height/2.;    // this flips the image at the same time
            float z = -height/(2.*tan(fov/2.));

            Vec3f dir = Vec3f(x, y, z);
            Vec3f orig = Vec3f(0,0,0);

            framebuffer[i+j*width] = cast_ray(orig, dir.normalize(), model, boxes, spheres, lights);
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

void renderGif(Model *model, const std::vector<Box> &boxes, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, const int &frame) {
    const int width    = 1024;
    const int height   = 768;
    const int fov      = M_PI/2.;
    std::vector<Vec3f> framebuffer(width*height);
    Vec3f eye = {0,0,0};
    Vec3f center = {0,0,-15};


    Matrix ModelView  = lookat(eye, center, Vec3f(0,1,0));
    Matrix Projection = Matrix::identity(4);
    Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);
    Projection[3][2] = -1.f/(eye-center).norm();

    Matrix z = (ViewPort*Projection*ModelView);

    //std::cout << ModelView << std::endl;

    float angle = 1 * nbFrame / 360.0;
    Matrix R = rot_z(angle * frame * M_PI / 180.0 );

    //std::cout << "apres rot" << std::endl;

    Matrix test = R * ModelView;
    
    //td::cout << test << std::endl;

    //std::cout << "rotation vecteur centre" << std::endl;

    Vec3f test2 = Vec3f(R * center);
    
    //std::cout << test2 << std::endl;

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            Vec3f center = {0,0,-15};
            //float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
            //float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            float x =  (i + 0.5) -  width/2.;
            float y = -(j + 0.5) + height/2.;    // this flips the image at the same time
            float z = -height/(2.*tan(fov/2.));

            Vec3f dir = Vec3f(x, y, z);
            Vec3f orig = Vec3f(0,0,0);

            dir = Vec3f(R*Matrix(dir));

            orig = center - Vec3f(R * center);

            //std::cout << orig << std::endl;


            framebuffer[i+j*width] = cast_ray(orig, dir.normalize(), model, boxes, spheres, lights);
        }
    }

    std::vector<unsigned char> pixmap(width*height*3);
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            pixmap[i*3+j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    std::string fileName = "tmp/out_" + std::to_string(frame) + ".jpg";
    stbi_write_jpg(fileName.c_str(), width, height, 3, pixmap.data(), 100);

    std::string frameString = "Gif frame " + std::to_string(frame) + "/" + std::to_string(nbFrame);
    std::cout << frameString << std::endl;
}

int main(int argc, char** argv) {
    int n = -1;
    unsigned char *pixmap = stbi_load("envmap.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    envmap = std::vector<Vec3f>(envmap_width*envmap_height);
    for (int j = envmap_height-1; j>=0 ; j--) {
        for (int i = 0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);

    //std::cout << envmap.size() << std::endl;

    Material      ivory(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3),   50.);
    Material      glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);

    Model *model = NULL;
    if (2==argc) {
        model = new Model(argv[1], glass);
    } else {
        model = new Model("duck.obj", glass);
    }

    std::vector<Box> boxes;
    boxes.push_back(Box(Vec3f(-5,    5,   -30), Vec3f(-2,   10,   -35), red_rubber));

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2,      glass));
    spheres.push_back(Sphere(Vec3f( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f( 7,    5,   -18), 4,     mirror));

    std::vector<Light>  lights;
    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));

    render(model, boxes, spheres, lights);

    
    system("mkdir -p tmp");
    for(size_t i = 0; i < nbFrame; i++){
        renderGif(model, boxes, spheres, lights, i);
    }
    system("ffmpeg -f image2 -framerate 24 -i tmp/out_%d.jpg out.gif");
    system("rm -r tmp");

    std::cout << "Gif créé" << std::endl;

    return 0;
}
