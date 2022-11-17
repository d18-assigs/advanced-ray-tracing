
o = newPlane(.05, .9, .35, .35, .55, .8, .75, 1, 1, 2);
Scale(o, 15, 15, 15);
RotateX(o, PI / 2);
Translate(o, 0, -4, 5);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newPlane(.05, .3, .35, .75, 1, 1, 1, 1, 1, 6);
Scale(o, (10/cos(PI/8)),11,11);
RotateY(o, PI / 8);
Translate(o, 10,-4, 10);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newPlane(.05, .3, .35, .75, 1, 1, 1, 1, 1, 6);
Scale(o, (10/cos(PI/8)),11,11);
RotateY(o, -PI / 8);
Translate(o, -10,-4, 10);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newSphere(.05, .95, .95, .13, 0,.18,.58, 1, 1, 6);
RotateX(o, PI / 2);
Translate(o, 0, 2, 7);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list); 
prev_object = o;
hierachycal_shpere(prev_object,2,.3,&object_list);

o = newCone(.05, .95, .35, .35, 1, .25, .25, 1, 1, 6);
Scale(o,1,1,5);
RotateX(o,-PI/2);
Translate(o,0,-4,7);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newCyl(.05, .95, .35, .35, 1, .25, .25, 1, 1, 6);
RotateZ(o, PI / 4);
RotateX(o, -PI / 2);
Translate(o,-2,-4,6.5);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newCyl(.05, .95, .35, .35, 1, .25, .25, 1, 1, 6);
RotateZ(o, PI / 4);
RotateX(o, -PI / 2);
Translate(o,2,-4,6.5);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

// p.px = 0;
// p.py = 25.5;
// p.pz = -3.5;
// p.pw = 1;
// l = newPLS(&p, .65, .65, .65);
// insertPLS(l, &light_list);

// p.px = -10;
// p.py = 10;
// p.pz = -25.5;
// p.pw = 1;
// l = newPLS(&p, .15, .15, .15);
// insertPLS(l, &light_list);

struct areaLS *aLS;
 o=newSphere(.05,.75,.05,.05,1,1,1,1,1,2);
 Scale(o,.8,.8,.8);
 Translate(o,0,6,1);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 aLS = newALS(o,500);
 insertALS(aLS,&aLS_list,&object_list);


