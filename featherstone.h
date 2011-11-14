#ifndef FEATHERSTONE_H
#define FEATHERSTONE_H

// featherstone

// FD
// ID (know loads)
// iterative inverse kinematics
// floating base
// closed systems
// constraints

vec FD( model &M, vec q, vec qd, vec tau);

mat bodypos( model &M, int b, vec q);

vec xformpos( mat X);

vec xformrot( mat X);

vec rel2abs( model &M, int b, vec q, vec relpos);

vec rel2abschainp( model &M, int b, vec q, vec relpos);

vec abs2rel( model &M, int b, vec q, vec abspos);

mat bodyjac( model &M, int b, vec q);

field<vec> IKpos( model &M, int b, mat Xd, vec q0);

#endif
