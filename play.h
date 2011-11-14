#ifndef PLAY_H
#define PLAY_H

#include <stdio.h>
#include <stdlib.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <armadillo>
#include "SDL.h"
//#include "obj.h"

/*
 * This code was created by Jeff Molofee '99 
 * (ported to Linux/SDL by Ti Leggett '01)
 *
 * If you've found this code useful, please let me know.
 *
 * Visit Jeff at http://nehe.gamedev.net/
 * 
 * or for port-specific comments, questions, bugreports etc. 
 * email to leggett@eecs.tulane.edu
 */
 
void Quit( int returnCode );

int resizeWindow( int width, int height );

void handleKeyPress( SDL_keysym *keysym );

int initGL();

int drawGLScene( model &M, field<vec> &pos, field<vec> &axis, field<vec> &endeff);

int play( model &M, spline &sp1, int saveImgs);

void prepDraw();

void drawCoord();

void drawBox( vec pos, vec axis, double w, double h, double d);

void drawBKGD();

void Screendump(char *destFile, short W, short H);

int play( model &M, field<vec> &Q, int saveImgs);

#endif
