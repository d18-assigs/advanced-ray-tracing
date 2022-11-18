 o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,6);		// Initialize a sphere
 Scale(o,1.5,.75,.75);					// Apply a few transforms (Translate * Rotate * Scale)
 RotateZ(o,PI/4);					
 Translate(o,2.0,2.5,1.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
 loadTexture(o,"./Texture/n_metal.ppm",2,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
 loadTexture(o,"./Texture/t_metal.ppm",1,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The 
  insertObject(o,&object_list);			// <-- If you don't insert the object into the object list,

 o=newSphere(.05,.95,.95,.75,0,.18,.58,1,1,6);
 Scale(o,.95,1.65,.65);
 RotateZ(o,-PI/1.5);
 Translate(o,-2.2,1.75,1.35);
 invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o,"./Texture/n_bricks.ppm",2,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
 loadTexture(o,"./Texture/t_bricks.ppm",1,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
 insertObject(o,&object_list);


 o=newSphere(.05,.95,.95,.75,0,.18,.58,1,1,6);
 Translate(o,0,0,2);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 loadTexture(o,"./Texture/a_geometric.pgm",3,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
loadTexture(o,"./Texture/n_fur.ppm",2,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
 loadTexture(o,"./Texture/t_fur.ppm",1,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
 insertObject(o,&object_list);

 o=newPlane(.05,.75,.05,.15,0.5,.29,.61,1,1,2);
 Scale(o,8,8,8);
 Translate(o,0,0,8);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

  o=newPlane(.05,.75,.05,.15,0.5,.29,.61,1,1,2);
 Scale(o,8,8,8);
 RotateX(o,-PI/2);
 Translate(o,0,8,0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

  o=newPlane(.05,.75,.05,.15,0.5,.29,.61,1,1,2);
 Scale(o,8,8,8);
 RotateX(o,PI/2);
 Translate(o,0,-8,0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 loadTexture(o,"./Texture/n_bricks_2.ppm",2,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
 loadTexture(o,"./Texture/t_wood_2.ppm",1,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
 insertObject(o,&object_list);

  o=newPlane(.05,.75,.05,.15,0.5,.29,.61,1,1,2);
 Scale(o,8,8,8);
 RotateY(o,-PI/2);
 Translate(o,-8,0,0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

  o=newPlane(.05,.75,.05,.15,0.5,.29,.61,1,1,2);
 Scale(o,8,8,8);
 RotateY(o,PI/2);
 Translate(o,8,0,0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);


 o=newSphere(.05,.75,.05,.05,1,1,1,1,1,2);
 Scale(o,.8,.8,.8);
 Translate(o,0,6,-5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 aLS = newALS(o,500);
 insertALS(aLS,&aLS_list,&object_list);