// o = newCyl(.05,.95,.95,.25,1,.18,.58,1,1.2,3);
// RotateX(o,-PI/2);
// Translate(o,0,1.55,3);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// insertObject(o,&object_list);
struct point3D up,cam_up;
up.px = 0;
up.py = 1;
up.pz = 0;
cam_up.px = 0;
cam_up.py = 0;
cam_up.pz = 1;
o = newPlane(0, 0 ,.95, .95,1,1,1, 1, 1, 2);
Scale(o, 8, 8, 8);
RotateX(o, PI / 2);
Translate(o, 0, -8, 0);
loadTexture(o, "Texture/coolScene/n_mirror.ppm", 2,
            &texture_list);  // 
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newSphere(.05, 0, .95, .25, 0, .18, .58, .2, 1.5, 3);
Scale(o,2,2,2);
Translate(o, -1, -7, 7);
invert(&o->T[0][0], &o->Tinv[0][0]);\
loadTexture(o,"Texture/coolScene/n_mirror.ppm",2,&texture_list);
// loadTexture(o,"./Texture/n_bricks.ppm",2,&texture_list);
// loadTexture(o,"./Texture/t_bricks.ppm",1,&texture_list);
insertObject(o, &object_list);

o = newSphere(.05, 0, .95, .25, .8, .18, .58, 1, 1.7, 3);
Translate(o, -1, -7, 4);
invert(&o->T[0][0], &o->Tinv[0][0]);
loadTexture(o,"Texture/coolScene/n_mirror.ppm",2,&texture_list);
// loadTexture(o,"./Texture/n_bricks.ppm",2,&texture_list);
// loadTexture(o,"./Texture/t_bricks.ppm",1,&texture_list);
insertObject(o, &object_list);


o = newSphere(.05, 0, .95, .95,0, .18, .58, 1, 1.1, 3);
Translate(o, 2, -7, 5.6);
invert(&o->T[0][0], &o->Tinv[0][0]);\
loadTexture(o,"Texture/n_fur.ppm",2,&texture_list);
// loadTexture(o,"./Texture/t_bricks.ppm",1,&texture_list);
insertObject(o, &object_list);

o = newCyl(.05, .05, .65, .65, .2,.2,.2, 1, 1, 9);
Scale(o,.7,.5,1);
RotateX(o,-PI/2);
loadTexture(o,"Texture/n_bricks.ppm",2,&texture_list);
loadTexture(o,"Texture/t_metal.ppm",1,&texture_list);
Translate(o, -4.3, -7, 5);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o,&object_list);


o = newCyl(.05, .05, .65, .65, .2,.2,.2, 1, 1, 9);
Scale(o,.7,.5,1);
RotateX(o,-PI/2);
loadTexture(o,"Texture/n_metal.ppm",2,&texture_list);
loadTexture(o,"Texture/t_metal.ppm",1,&texture_list);
Translate(o, 4.3, -7, 5);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o,&object_list);


o = newPlane(.05, .75,.1,.1, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
Translate(o, 0, 0, 8);
invert(&o->T[0][0], &o->Tinv[0][0]);
RotateZ(o,PI);
loadTexture(o, "Texture/t_robosoccer.ppm", 1,
            &texture_list); 
            loadTexture(o, "Texture/coolScene/n_paper.ppm", 2,
            &texture_list);  // This 
insertObject(o, &object_list);

o = newPlane(.05, .75, 0,0, .6,.6,.6, 1, 1, 1);
Scale(o, 8, 8, 8);
RotateX(o, -PI / 2);
Translate(o, 0, 8, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
            loadTexture(o, "Texture/coolScene/n_paper.ppm", 2,
            &texture_list);  // This 
insertObject(o, &object_list);



o = newPlane(.05, .25, .15, .15, .8,.8,.8, 1, 1, 1);
Scale(o, 8, 8, 8);
RotateY(o, -PI / 2);
Translate(o, -8, 0, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
loadTexture(o, "Texture/coolScene/belly.ppm", 1,
            &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
loadTexture(o, "Texture/n_wood.ppm", 2,
            &texture_list); 
// loadTexture(o, "Texture/t_bricks.ppm", 1,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .05, .8,.8,.8, 1, 1, 1);
Scale(o, 8, 8, 8);
RotateY(o, PI / 2);
RotateX(o,PI);
Translate(o, 8, 0, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
loadTexture(o, "Texture/coolScene/top.ppm", 1,
            &texture_list);  // This 
loadTexture(o, "Texture/n_wood.ppm", 2,
            &texture_list); 
insertObject(o, &object_list);




// o = newPlane(.05, .75, .05, .05, 1, 1, 1, 1, 1, 2);
// up.px = 0;
// up.py = 0;
// up.pz = -1;
// cam_up.px = 1;
// cam_up.py = -1;
// cam_up.pz = -1;
// RotateObjToVec(o,&up,&cam_up);
// Translate(o, -4, 6, 5);
// invert(&o->T[0][0], &o->Tinv[0][0]);
// aLS = newALS(o, 500);
// insertALS(aLS, &aLS_list, &object_list);



o = newCyl(.05, .75, .05, .05,1,1,1, 1, 1, 2);
// o = newCyl(.05, .75, .05, .05, .92, 0, .48, 1, 1, 2);
up.px = 0;
up.py = 1;
up.pz = 0;
cam_up.px = 0;
cam_up.py = 0;
cam_up.pz = 1;
Scale(o,.3,.3,3);
RotateObjToVec(o,&up,&cam_up);
Translate(o, 4.3, -5, 5);
invert(&o->T[0][0], &o->Tinv[0][0]);
aLS = newALS(o, 700);
insertALS(aLS, &aLS_list, &object_list);




o = newCyl(.05, .75, .05, .05, 1,1,1, 1, 1, 2);
// o = newCyl(.05, .75, .05, .05, .43, .12, .52, 1, 1, 2);
up.px = 0;
up.py = 1;
up.pz = 0;
cam_up.px = 0;
cam_up.py = 0;
cam_up.pz = 1;
Scale(o,.3,.3,3);
RotateObjToVec(o,&up,&cam_up);
Translate(o, -4.3, -5, 5);
invert(&o->T[0][0], &o->Tinv[0][0]);
aLS = newALS(o, 700);
insertALS(aLS, &aLS_list, &object_list);
