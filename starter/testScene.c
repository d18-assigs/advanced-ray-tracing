// o = newCyl(.05,.95,.95,.25,1,.18,.58,1,1.2,3);
// RotateX(o,-PI/2);
// Translate(o,0,1.55,3);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// insertObject(o,&object_list);

o = newSphere(.05, .95, .95, .25, 0, .18, .58, .4, 1.2, 3);
Translate(o, 0, 1.55, 1);
invert(&o->T[0][0], &o->Tinv[0][0]);
loadTexture(o,"./Texture/n_bricks.ppm",2,&texture_list);
loadTexture(o,"./Texture/t_bricks.ppm",1,&texture_list);
loadTexture(o,"./Texture/a_geometric.pgm",1,&texture_list);
insertObject(o, &object_list);
prev_object = o;
hierachycal_shpere(prev_object, 2, .3, &object_list, &texture_list);

o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
Translate(o, 0, 0, 8);
invert(&o->T[0][0], &o->Tinv[0][0]);
// loadTexture(o, "./Texture/n_bricks.ppm", 2,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// loadTexture(o, "./Texture/t_bricks.ppm", 1,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
RotateX(o, -PI / 2);
Translate(o, 0, 8, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
// loadTexture(o, "./Texture/n_bricks.ppm", 2,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// loadTexture(o, "./Texture/t_bricks.ppm", 1,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
RotateX(o, PI / 2);
Translate(o, 0, -8, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
loadTexture(o, "./Texture/n_metal.ppm", 2,
            &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
loadTexture(o, "./Texture/t_metal.ppm", 1,
            &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);

// Robosoccer
o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
// Scale(o, 8, 8, 8);
RotateX(o, PI / 2);
Translate(o, 0, -3, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
loadTexture(o, "./Texture/t_robosoccer.ppm", 1,
            &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
RotateY(o, -PI / 2);
Translate(o, -8, 0, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
// loadTexture(o, "./Texture/n_bricks.ppm", 2,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// loadTexture(o, "./Texture/t_bricks.ppm", 1,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .15, 0.5, .29, .61, 1, 1, 2);
Scale(o, 8, 8, 8);
RotateY(o, PI / 2);
Translate(o, 8, 0, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
// loadTexture(o, "./Texture/n_bricks.ppm", 2,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
// loadTexture(o, "./Texture/t_bricks.ppm", 1,
//             &texture_list);  // This loads a texture called 'mosaic2.ppm'. The
insertObject(o, &object_list);


o = newSphere(.05, .75, .05, .05, 1, 1, 1, 1, 1, 2);
Scale(o, .8, .8, .8);
Translate(o, 0, 6, -5);
invert(&o->T[0][0], &o->Tinv[0][0]);
aLS = newALS(o, 500);
insertALS(aLS, &aLS_list, &object_list);