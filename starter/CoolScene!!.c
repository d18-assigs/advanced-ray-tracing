// o = newCyl(.05,.95,.95,.25,1,.18,.58,1,1.2,3);
// RotateX(o,-PI/2);
// Translate(o,0,1.55,3);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// insertObject(o,&object_list);

o = newSphere(.05, 0, .95, .25, 0, .18, .58, .4, 1.2, 3);
Translate(o, 0, 0, 1);
invert(&o->T[0][0], &o->Tinv[0][0]);\
// loadTexture(o,"./Texture/n_bricks.ppm",2,&texture_list);
// loadTexture(o,"./Texture/t_bricks.ppm",1,&texture_list);
insertObject(o, &object_list);
// prev_object = o;
// fractal_sphere(prev_object, 2, .3, &object_list, &texture_list);

o = newPlane(.05, .75,0,0, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
Translate(o, 0, 0, 8);
invert(&o->T[0][0], &o->Tinv[0][0]);
RotateZ(o,PI);
// loadTexture(o, "Texture/t_robosoccer.ppm", 1,
//             &texture_list); 
loadTexture(o, "Texture/coolScene/n_paper.ppm", 2,
            &texture_list); 
insertObject(o, &object_list);

// o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
// Scale(o, 8, 8, 8);
// RotateX(o, -PI / 2);
// Translate(o, 0, 8, 0);
// invert(&o->T[0][0], &o->Tinv[0][0]);
// // loadTexture(o, "Texture/n_bricks.ppm", 2,
// //             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// // loadTexture(o, "Texture/t_bricks.ppm", 1,
// //             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// insertObject(o, &object_list);

// o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
// Scale(o, 8, 8, 8);
// RotateX(o, PI / 2);
// Translate(o, 0, -8, 0);
// invert(&o->T[0][0], &o->Tinv[0][0]);
// // loadTexture(o, "./Texture/n_metal.ppm", 2,
// //             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// // loadTexture(o, "./Texture/t_metal.ppm", 1,
// //             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// insertObject(o, &object_list);

// Robosoccer
o = newPlane(0, 0 ,.95, .95, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
RotateX(o, PI / 2);
Translate(o, 0, -8, 0);
loadTexture(o, "./Texture/coolScene/n_mirror.ppm", 2,
            &texture_list);  // 
invert(&o->T[0][0], &o->Tinv[0][0]);

 // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
RotateY(o, -PI / 2);
Translate(o, -8, 0, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
// loadTexture(o, "Texture/n_bricks.ppm", 2,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// loadTexture(o, "Texture/t_bricks.ppm", 1,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
RotateY(o, PI / 2);
Translate(o, 8, 0, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
// loadTexture(o, "Texture/n_bricks.ppm", 2,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// loadTexture(o, "Texture/t_bricks.ppm", 1,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);

struct point3D up,point_to;
//  p.px=0;
//  p.py=0;
//  p.pz=-20;
//  p.pw=1;
//  l=newPLS(&p,.5,.5,.5);
//  insertPLS(l,&light_list);
// o = newPlane(.05, .75, .05, .05, 1, 1, 1, 1, 1, 2);
// up.px = 0;
// up.py = 0;
// up.pz = -1;
// point_to.px = 1;
// point_to.py = -1;
// point_to.pz = -1;
// RotateObjToVec(o,&up,&point_to);
// Translate(o, -4, 6, 5);
// invert(&o->T[0][0], &o->Tinv[0][0]);
// aLS = newALS(o, 500);
// insertALS(aLS, &aLS_list, &object_list);

// o = newPlane(.05, .75, .05, .05, 1, 1, 1, 1, 1, 2);
// up.px = 0;
// up.py = 0;
// up.pz = 1;
// point_to.px = 1;
// point_to.py = -1;
// point_to.pz = -1;
// RotateObjToVec(o,&up,&point_to);
// Translate(o, -4, 6, 5);
// invert(&o->T[0][0], &o->Tinv[0][0]);
// aLS = newALS(o, 500);
// insertALS(aLS, &aLS_list, &object_list);


o = newCyl(.05, .75, .05, .05,1,1,1, 1, 1, 2);
// o = newCyl(.05, .75, .05, .05, .92, 0, .48, 1, 1, 2);
up.px = 0;
up.py = 1;
up.pz = 0;
point_to.px = 0;
point_to.py = 0;
point_to.pz = 1;
Scale(o,.2,.2,3);

RotateObjToVec(o,&up,&point_to);
Translate(o, 4, -5, 4);
invert(&o->T[0][0], &o->Tinv[0][0]);
aLS = newALS(o, 500);
insertALS(aLS, &aLS_list, &object_list);




o = newCyl(.05, .75, .05, .05, 1,1,1, 1, 1, 2);
// o = newCyl(.05, .75, .05, .05, .43, .12, .52, 1, 1, 2);
up.px = 0;
up.py = 1;
up.pz = 0;
point_to.px = 0;
point_to.py = 0;
point_to.pz = 1;
Scale(o,.2,.2,3);
RotateObjToVec(o,&up,&point_to);
Translate(o, -4, -5, 4);
invert(&o->T[0][0], &o->Tinv[0][0]);
aLS = newALS(o, 500);
insertALS(aLS, &aLS_list, &object_list);