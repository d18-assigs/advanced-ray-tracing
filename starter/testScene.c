

o = newCyl(.05,.95,.95,.25,1,.18,.58,1,1.2,3);
RotateX(o,-PI/2);
Translate(o,0,1.55,3);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

o=newSphere(.05,.95,.95,.25,0,.18,.58,.4,1.2,3);
Translate(o,0,1.55,1);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);
prev_object = o;
hierachycal_shpere(prev_object,2,.3,&object_list);



//   o=newPlane(.05,.75,.05,.15,0.5,.29,.61,1,1,2);
//  Scale(o,22,22,22);
//  RotateX(o,PI/2);
//  Translate(o,0,-8,22);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);


// //  o=newSphere(.05,.95,.95,.75,0,.18,.58,1,1,6);
// //  Translate(o,0,0,2);
// //  invert(&o->T[0][0],&o->Tinv[0][0]);
// //  insertObject(o,&object_list);

 o=newPlane(.05,.75,.05,.15,0.5,.29,.61,1,1,2);
 Scale(o,8,8,8);
 Translate(o,0,0,8);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);


//   o=newPlane(.05,.75,.05,.05,0.5,.29,.61,1,1,2);
//  Scale(o,8,8,8);
//  RotateY(o,-PI/2);
//  Translate(o,-8,0,0);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//   o=newPlane(.05,.75,.05,.05,0.5,.29,.61,1,1,2);
//  Scale(o,8,8,8);
//  RotateY(o,PI/2);
//  Translate(o,8,0,0);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//   o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
//  Scale(o,11,11,11);
//  RotateZ(o,PI/4);
//  RotateX(o,PI/2);
//  Translate(o,0,-4,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

 // Insert a single point light source. We set up its position as a point structure, and specify its
 // colour in terms of RGB (in [0,1]).

 p.px=0;
 p.py=20;
 p.pz=-5;
 p.pw=1;
 l=newPLS(&p,1,1,1);
 insertPLS(l,&light_list);


//  o=newSphere(.05,.75,.05,.05,1,1,1,1,1,2);
//  Scale(o,.8,.8,.8);
//  Translate(o,0,6,6);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  aLS = newALS(o,500);
//  insertALS(aLS,&aLS_list,&object_list);

//   o = newCyl(.05,.75,.05,.05,1,1,1,1,1,2);
//   Translate(o,0,0,-.5);
//   Scale(o,.3,.3,3);
//   RotateX(o,-PI/3);
//   Translate(o,0,1,2);
//   invert(&o->T[0][0],&o->Tinv[0][0]);
//   aLS = newALS(o,1000);
//   insertALS(aLS,&aLS_list,&object_list);


  // o = newCyl(.05,.75,.05,.05,1,1,1,1,1,2);
  // Translate(o,0,0,-.5);
  // Scale(o,.1,.1,5);
  // // RotateY(o,-PI/6);
  // RotateX(o,PI/2);
  // Translate(o,5.5,0,3);
  // invert(&o->T[0][0],&o->Tinv[0][0]);
  // aLS = newALS(o,100);
  // insertALS(aLS,&aLS_list,&object_list);


  // o = newCyl(.05,.75,.05,.05,1,1,1,1,1,2);
  // Translate(o,0,0,-.5);
  // Scale(o,.1,.1,2);
  // RotateY(o,-PI/6);
  // // RotateX(o,PI/2);
  // Translate(o,0,2,2);
  // invert(&o->T[0][0],&o->Tinv[0][0]);
  // aLS = newALS(o,500);
  // insertALS(aLS,&aLS_list,&object_list);
  
  
  // o=newSphere(.05,.75,.05,.05,1,1,1,1,1,2);
  // RotateX(o,-PI/2);
  // Translate(o,0,9,7);
  // invert(&o->T[0][0],&o->Tinv[0][0]);
  // aLS = newALS(o,500);
  // insertALS(aLS,&aLS_list,&object_list);

 // End of simple scene for Assignment 2
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
 
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
 // TO DO: For Assignment 3 you *MUST* define your own cool scene.
 //	   We will be looking for the quality of your scene setup, the use of hierarchical or composite
 //	   objects that are more interesting than the simple primitives from A2, the use of textures
 //        and other maps, illumination and illumination effects such as soft shadows, reflections and
 //        transparency, and the overall visual quality of your result. Put some work into thinking
 //        about these elements when designing your scene.
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
