#include <stdio.h>
#include <stdlib.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <armadillo>
#include <algorithm>
#include "SDL.h"
#include <IL/il.h>
#include "rb.h"
#include "rbbox.h"
#include "rbobj.h"
#include "model.h"
#include "featherstone.h"
#include "spline.h"

#include "play.h"

using namespace std;
using namespace arma;

/* screen width, height, and bit depth */
#define SCREEN_WIDTH  720
#define SCREEN_HEIGHT 540
#define SCREEN_BPP     16

/* Define our booleans */
#define TRUE  1
#define FALSE 0

/* This is our SDL surface */
SDL_Surface *surface;

GLUquadricObj *quadratic;

int done;

GLint tidx = 0;
int maxTidx = 1;
GLfloat sec = 1;
GLfloat scenex = 0;
GLfloat sceney = 0;
GLfloat scenez = 0;
GLfloat offsetz = 
1.4;

GLfloat sceneroty = 0;
GLfloat scenerotx = 0;

wall_clock timer;
GLfloat durationActive;
GLfloat durationInactive;

char imgname[50];


/* Here goes our drawing code */
int drawGLScene( model &M, field<vec> &pos, field<vec> &axis, field<vec> &endeff)
{
    /* rotational vars for the triangle and quad, respectively */
    static GLfloat rtri, rquad, time;
    /* These are to calculate our fps */
    static GLint T0     = 0;
    static GLint Frames = 0;

	int N = M.getN();
	
	prepDraw();
	
	drawCoord();
	
	for (int i = 0; i < N; i++)
	{
		prepDraw();
		
		M.getRB(i)->draw( pos(tidx,i), axis(tidx,i));
		
		drawCoord();
	}
	
	prepDraw();
	
	glColor3f( 0.0f, 0.0f, 1);
	glBegin(GL_LINE_STRIP);
	for (int i = 1; i <= tidx; i++)
	{
		glVertex3f( endeff(i)(0), endeff(i)(1), endeff(i)(2) );
	}
	glEnd();

    prepDraw();

//	drawBKGD();

    // CHANGE THIS
//	scenex = -pos(tidx,N-1)(0);
   sceney = -pos(tidx,N-1)(1);
//    scenez = -pos(tidx,N-1)(2) - offsetz*M.getN();
	sceneroty = -double(tidx)*180/maxTidx;

    
    /* Draw it to the screen */
    SDL_GL_SwapBuffers( );
    
    /* Gather our frames per second */
    Frames++;
    {
		GLint t = SDL_GetTicks();
		if (t - T0 >= 5000)
		{
			GLfloat seconds = (t - T0) / 1000.0;
			GLfloat fps = Frames / seconds;
			printf("%d frames in %g seconds = %g FPS\n", Frames, seconds, fps);
			T0 = t;
			Frames = 0;
		}
    }
    

	rtri+= 0.2f;
	rquad+= 0.15f;

    return( TRUE );
}

/* function to release/destroy our resources and restoring the old desktop */
void Quit( int returnCode )
{
	done = FALSE;
	
	/* delete quadric */
	gluDeleteQuadric( quadratic );
    /* clean up the window */    
    SDL_Quit( );

    /* and exit appropriately */
    exit( returnCode );
}

/* function to reset our viewport after a window resize */
int resizeWindow( int width, int height )
{
    /* Height / width ration */
    GLfloat ratio;
 
    /* Protect against a divide by zero */
   if ( height == 0 )
	height = 1;

    ratio = ( GLfloat )width / ( GLfloat )height;

    /* Setup our viewport. */
    glViewport( 0, 0, ( GLsizei )width, ( GLsizei )height );

    /* change to the projection matrix and set our viewing volume. */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity(); // reset

    // viewing angle, aspect ratio, closest obj, farthest viewing depth
    gluPerspective(45.0f,ratio,0.1f,200.0f);

    // modelview holds object information
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
	glTranslatef(0.0f,0.0f,-25.0f);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glColor3f(1.0f,1.0f,1.0f);

    return( TRUE );
}

/* function to handle key press events */
void handleKeyPress( SDL_keysym *keysym )
{
    switch ( keysym->sym )
	{
	case SDLK_ESCAPE:
	    /* ESC key was pressed */
	    Quit( 0 );
	    break;
	case SDLK_PAGEUP:
		scenez += 0.2f;
		break;
	case SDLK_PAGEDOWN:
		scenez -= 0.2f;
		break;
	case SDLK_w:
		sceney += 0.2f;
		break;
	case SDLK_s:
		sceney -= 0.2f;
		break;
	case SDLK_d:
		scenex += 0.2f;
		break;
	case SDLK_a:
		scenex -= 0.2f;
		break;
	case SDLK_RIGHT:
		sceneroty += 5.0f;
		break;
	case SDLK_LEFT:
		sceneroty -= 5.0f;
		break;
	case SDLK_UP:
		scenerotx += 5.0f;
		break;
	case SDLK_DOWN:
		scenerotx -= 5.0f;
		break;
	case SDLK_r:
		tidx = 0;
		sec = 1;
		timer.tic();
		durationInactive = 0.0;
		break;
	case SDLK_i:
		tidx = min(tidx+1,maxTidx-1);
		cout << "iteration: " << tidx << endl;
		break;
	case SDLK_u:
		tidx = max(tidx-1,0);
		cout << "iteration: " << tidx << endl;
		break;
	case SDLK_p:
		sprintf(imgname,"animpi.tga");
		Screendump( imgname, SCREEN_WIDTH, SCREEN_HEIGHT);
			   			    
	   	ilLoadImage(imgname);
	   	sprintf(imgname,"animpi.png");
	   	ilSave(IL_PNG,imgname);
	   	break;
	case SDLK_F1:
	    /* F1 key was pressed
	     * this toggles fullscreen mode
	     */
	    SDL_WM_ToggleFullScreen( surface );
	    break;
	default:
	    break;
	}

    return;
}

/* general OpenGL initialization function */
int initGL()
{

    /* Enable smooth shading */
    glShadeModel( GL_SMOOTH );

    /* Set the background black */
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );

    /* Depth buffer setup */
    glClearDepth( 1.0f );

    /* Enables Depth Testing */
    glEnable( GL_DEPTH_TEST );
    
    /* Enable alpha */
//    glEnable( GL_LIGHTING );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    /* The Type Of Depth Test To Do */
    glDepthFunc( GL_LEQUAL );

    /* Really Nice Perspective Calculations */
    glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );
    
    /* quadric object */
    quadratic = gluNewQuadric();    
    
//    gluQuadricNormals( quadratic, GLU_SMOOTH );
// gluQuadricTexture( quadratic, GL_TRUE );

    return( TRUE );
}

int play( model &M, spline &sp1, int saveImgs)
{
    /* Flags to pass to SDL_SetVideoMode */
    int videoFlags;
    /* main loop variable */
    done = FALSE;
    /* used to collect events */
    SDL_Event event;
    /* this holds some info about our display */
    const SDL_VideoInfo *videoInfo;
    /* whether or not the window is active */
    int isActive = TRUE;
    /* whether or not the window was active last loop */
    int wasAcive = FALSE;
    /* save time at which inactivity started */
    GLfloat timeAtInactive = 0.0;
    /* updated once a loop, holds the time (s) active */
    durationActive = 0.0;
    /* duration to discard */
    durationInactive = 0.0;

    /* initialize SDL */
    if ( SDL_Init( SDL_INIT_VIDEO ) < 0 )
	{
	    fprintf( stderr, "Video initialization failed: %s\n",
		     SDL_GetError( ) );
	    Quit( 1 );
	}

    /* Fetch the video info */
    videoInfo = SDL_GetVideoInfo( );

    if ( !videoInfo )
	{
	    fprintf( stderr, "Video query failed: %s\n",
		     SDL_GetError( ) );
	    Quit( 1 );
	}

    /* the flags to pass to SDL_SetVideoMode */
    videoFlags  = SDL_OPENGL;          /* Enable OpenGL in SDL */
    videoFlags |= SDL_GL_DOUBLEBUFFER; /* Enable double buffering */
    videoFlags |= SDL_HWPALETTE;       /* Store the palette in hardware */
    videoFlags |= SDL_RESIZABLE;       /* Enable window resizing */

    /* This checks to see if surfaces can be stored in memory */
    if ( videoInfo->hw_available )
	videoFlags |= SDL_HWSURFACE;
    else
	videoFlags |= SDL_SWSURFACE;

    /* This checks if hardware blits can be done */
    if ( videoInfo->blit_hw )
	videoFlags |= SDL_HWACCEL;

    /* Sets up OpenGL double buffering */
    SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );

    /* get a SDL surface */
    surface = SDL_SetVideoMode( SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP,
				videoFlags );

    /* Verify there is a surface */
    if ( !surface )
	{
	    fprintf( stderr,  "Video mode set failed: %s\n", SDL_GetError( ) );
	    Quit( 1 );
	}

    /* initialize OpenGL */
    initGL( );

    /* resize the initial window */
    resizeWindow( SCREEN_WIDTH, SCREEN_HEIGHT );

	/* start timer */
	timer.tic();
	
	scenez = -offsetz*M.getN();

	// CALCULATE POS AND AXIS
	mat q = M.getQ();
	int nsteps = q.n_rows;
	mat X;
	field<vec> pos( nsteps , M.getN() );
	field<vec> axis( nsteps , M.getN() );
	field<vec> endeff( nsteps );
	field<vec> endeffk( nsteps );
	vec RR = zeros<vec>(3);
	RR(0) = 1.0;

	for (int i = 0; i < nsteps; i++)
	{
		for (int  j = 0; j < M.getN(); j++)
		{
//			X = bodypos( M, j, sp1.eval( M.getT()(i) ) );
			X = bodypos( M, j, vec( trans( q.row(i) ) ) );
			pos(i,j) = xformpos( X);
			axis(i,j) = xformrot( X);
		}	
		endeff(i) = rel2abs( M, M.getN()-1, vec(trans(q.row(i))), RR );
		endeffk(i) = rel2abs( M, M.getN()-1, sp1.eval( M.getT()(i) ), RR );
	}
	
	maxTidx = M.getT().n_elem-1;
	
	/* initialize IL */
	ilInit();
	
    /* wait for events */ 
    while ( !done )
	{
	    /* handle the events in the queue */

	    while ( SDL_PollEvent( &event ) )
		{
		    switch( event.type )
			{
			case SDL_ACTIVEEVENT:
			    /* Something's happend with our focus
			     * If we lost focus or we are iconified, we
			     * shouldn't draw the screen
			     */
			    if ( event.active.gain == 0 )
				{
					isActive = FALSE;
					timeAtInactive = timer.toc();
				}
			    else
				{
					isActive = TRUE;
					durationInactive += timer.toc() - timeAtInactive;
				}
			    break;			    
			case SDL_VIDEORESIZE:
			    /* handle resize event */
			    surface = SDL_SetVideoMode( event.resize.w,
							event.resize.h,
							16, videoFlags );
			    if ( !surface )
				{
				    fprintf( stderr, "Could not get a surface after resize: %s\n", SDL_GetError( ) );
				    Quit( 1 );
				}
			    resizeWindow( event.resize.w, event.resize.h );
			    break;
			case SDL_KEYDOWN:
			    /* handle key presses */
			    handleKeyPress( &event.key.keysym );
			    break;
			case SDL_QUIT:
			    /* handle quit requests */
			    done = TRUE;
			    break;
			default:
			    break;
			}
		}

		durationActive = timer.toc() - durationInactive;

	    /* draw the scene */
	    if ( isActive )
		{
			if ( durationActive < M.getT()(tidx) ){ /* wait */}
			else
			{
				if ( durationActive > sec )
				{
					cout << "elapsed time: " << sec << " s" << endl;
					sec++;
				}
				
			    /* Clear The Screen And The Depth Buffer */
			    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

				prepDraw();
				glColor3f( 1, 0.0f, 0.0f);
				glBegin(GL_LINE_STRIP);
				for (int i = 1; i <= maxTidx; i++)
				{
					glVertex3f( endeffk(i)(0), endeffk(i)(1), endeffk(i)(2) );
				}
				glEnd();
				glColor3f( 0.0f, 1, 0.0f);
				glBegin(GL_LINES);
					glVertex3f( endeff(tidx)(0), endeff(tidx)(1), endeff(tidx)(2) );
					glVertex3f( endeffk(tidx)(0), endeffk(tidx)(1), endeffk(tidx)(2) );
				glEnd();

				// if time is greater than sp1.getT()(countr)
				// drawLinks using GL_LINE_STRIP
				drawGLScene( M, pos, axis, endeff); // YAYYYYYYYY

				if ( saveImgs && tidx != maxTidx )
				{
					sprintf(imgname,"anim%d.tga",tidx);
					Screendump( imgname, SCREEN_WIDTH, SCREEN_HEIGHT);
				   			    
				   	ilLoadImage(imgname);
				   	sprintf(imgname,"anim%d.png",tidx);
				   	ilSave(IL_PNG,imgname);
			   	}
			   	
				tidx = min((int)tidx+1,maxTidx);
			}
		}
	}

    /* clean ourselves up and exit */
    Quit( 0 );

    /* Should never get here */
    return( 0 );
}

void prepDraw()
{
	glLoadIdentity();
	glTranslatef(scenex,sceney,scenez);
	glRotatef(sceneroty,0.0f,1.0f,0.0f);
	glRotatef(-sceneroty,0.0f,1.0f,0.0f);
	glRotatef(scenerotx,1.0f,0.0f,0.0f);
	glRotatef(sceneroty,0.0f,1.0f,0.0f);
}

void drawCoord()
{	
	glColor3f(0.7f,0.0f,0.0f);
	glRotatef( 90, 0.0f, 1, 0.0f);
	gluCylinder( quadratic , 0.025, 0.025, 0.5, 32 , 32);
	glTranslatef( 0.0f, 0.0f, 0.5);
	gluCylinder( quadratic , 0.05, 0.0, 0.10, 32 , 32);
	glTranslatef( 0.0f, 0.0f, -0.5);
	//	gluSphere( quadratic , 0.25f, 32 , 32 );
	
	glColor3f(0.0f,0.7f,0.0f);
	glRotatef( -90, 1, 0.0f, 0.0f);
	gluCylinder( quadratic , 0.025, 0.025, 0.5, 32 , 32);
	glTranslatef( 0.0f, 0.0f, 0.5);
	gluCylinder( quadratic , 0.05, 0.0, 0.10, 32 , 32);
	glTranslatef( 0.0f, 0.0f, -0.5);
	
	glColor3f(0.0f,0.0f,0.7f);
	glRotatef( -90, 0.0f, 1, 0.0f);
	gluCylinder( quadratic , 0.025, 0.025, 0.5, 32 , 32);
	glTranslatef( 0.0f, 0.0f, 0.5);
	gluCylinder( quadratic , 0.05, 0.0, 0.10, 32 , 32);
	glTranslatef( 0.0f, 0.0f, -0.5);

}

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
 
void drawBKGD()
{
//    glTranslatef(0.0,0.0,offsetz/2*M.getN());
    glColor4f(1.0,0.0,0.0,0.3f);
    
    float L = 5;
    glBegin(GL_QUADS);
		glVertex3f(0, -L, -L);
		glVertex3f(0,  L, -L);
		glVertex3f(0,  L,  L);
		glVertex3f(0, -L, L);
    glEnd();

	glColor4f(0.0,1.0,0.0,0.3f);
    glBegin(GL_QUADS);
		glVertex3f(0, 0, -L);
		glVertex3f( L, 0, -L);
		glVertex3f( L, 0,  L);
		glVertex3f(0, 0,  L);
    glEnd();

	glColor4f(0.0,0.0,1.0,0.3f);
    glBegin(GL_QUADS);
		glVertex3f(0, -L, 0);
		glVertex3f( L, -L, 0);
		glVertex3f( L, L,  0);
		glVertex3f(0, L, 0);
    glEnd();
    
    
    glColor4f(0.0,0.6,1.0,0.5f);

    float w = 3;
    float h = 4;
    float d = 7; 

   	glTranslatef(2.8f + w, 0, -d/2);
    // face 1
    glBegin(GL_QUADS);
		glVertex3f(-w, -h, 0);
		glVertex3f(-w,  h, 0);
		glVertex3f(-w,  h, d);
		glVertex3f(-w, -h, d);
    glEnd();
    
    // face 2    
    glBegin(GL_QUADS);
		glVertex3f( w, -h, 0);
		glVertex3f( w,  h, 0);
		glVertex3f( w,  h, d);
		glVertex3f( w, -h, d);
    glEnd();
    
    // face 3    
    glBegin(GL_QUADS);
		glVertex3f(-w, -h, 0);
		glVertex3f(-w, -h, d);
		glVertex3f( w, -h, d);
		glVertex3f( w, -h, 0);
    glEnd();

	// face 4
    glBegin(GL_QUADS);
		glVertex3f(-w,  h, 0);
		glVertex3f(-w,  h, d);
		glVertex3f( w,  h, d);
		glVertex3f( w,  h, 0);
    glEnd();
    
    // face 5
    glBegin(GL_QUADS);
		glVertex3f(-w, -h, d);
		glVertex3f(-w,  h, d);
		glVertex3f( w,  h, d);
		glVertex3f( w, -h, d);
    glEnd();
    
    // face 6
    glBegin(GL_QUADS);
		glVertex3f(-w, -h, 0);
		glVertex3f(-w,  h, 0);
		glVertex3f( w,  h, 0);
		glVertex3f( w, -h, 0);
    glEnd();

}
void Screendump(char *destFile, short W, short H)
{
	FILE   *out = fopen(destFile, "w");
	char   pixel_data[3*W*H];
	short  TGAhead[] = {0, 2, 0, 0, 0, 0, W, H, 24};
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, W, H, GL_BGR, GL_UNSIGNED_BYTE, pixel_data);
	fwrite(&TGAhead, sizeof(TGAhead), 1, out);
	fwrite(pixel_data, 3*W*H, 1, out); fclose(out);
}/*       ______________________________________________________________   creates a 24-bit uncompressed true color tga-file, width W and height   H. for a complete description of tga go to www.wotsit.org. if you use   a doublebuffered configuration, be sure to swap buffers before you   call the screendump function.*/


































int play( model &M, field<vec> &Q, int saveImgs)
{
    /* Flags to pass to SDL_SetVideoMode */
    int videoFlags;
    /* main loop variable */
    done = FALSE;
    /* used to collect events */
    SDL_Event event;
    /* this holds some info about our display */
    const SDL_VideoInfo *videoInfo;
    /* whether or not the window is active */
    int isActive = TRUE;
    /* whether or not the window was active last loop */
    int wasAcive = FALSE;
    /* save time at which inactivity started */
    GLfloat timeAtInactive = 0.0;
    /* updated once a loop, holds the time (s) active */
    durationActive = 0.0;
    /* duration to discard */
    durationInactive = 0.0;

    /* initialize SDL */
    if ( SDL_Init( SDL_INIT_VIDEO ) < 0 )
	{
	    fprintf( stderr, "Video initialization failed: %s\n",
		     SDL_GetError( ) );
	    Quit( 1 );
	}

    /* Fetch the video info */
    videoInfo = SDL_GetVideoInfo( );

    if ( !videoInfo )
	{
	    fprintf( stderr, "Video query failed: %s\n",
		     SDL_GetError( ) );
	    Quit( 1 );
	}

    /* the flags to pass to SDL_SetVideoMode */
    videoFlags  = SDL_OPENGL;          /* Enable OpenGL in SDL */
    videoFlags |= SDL_GL_DOUBLEBUFFER; /* Enable double buffering */
    videoFlags |= SDL_HWPALETTE;       /* Store the palette in hardware */
    videoFlags |= SDL_RESIZABLE;       /* Enable window resizing */

    /* This checks to see if surfaces can be stored in memory */
    if ( videoInfo->hw_available )
	videoFlags |= SDL_HWSURFACE;
    else
	videoFlags |= SDL_SWSURFACE;

    /* This checks if hardware blits can be done */
    if ( videoInfo->blit_hw )
	videoFlags |= SDL_HWACCEL;

    /* Sets up OpenGL double buffering */
    SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );

    /* get a SDL surface */
    surface = SDL_SetVideoMode( SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP,
				videoFlags );

    /* Verify there is a surface */
    if ( !surface )
	{
	    fprintf( stderr,  "Video mode set failed: %s\n", SDL_GetError( ) );
	    Quit( 1 );
	}

    /* initialize OpenGL */
    initGL( );

    /* resize the initial window */
    resizeWindow( SCREEN_WIDTH, SCREEN_HEIGHT );

	/* start timer */
	timer.tic();
	
	scenez = -offsetz*M.getN();

	// CALCULATE POS AND AXIS
	maxTidx = Q.n_elem;
	mat X;
	field<vec> pos( maxTidx , M.getN() );
	field<vec> axis( maxTidx , M.getN() );
	field<vec> endeff( maxTidx );
	vec RR = zeros<vec>(3);
	RR(0) = 1.0;

	for (int i = 0; i < maxTidx; i++)
	{
		for (int  j = 0; j < M.getN(); j++)
		{
			X = bodypos( M, j, Q(i) );
			pos(i,j) = xformpos( X);
			axis(i,j) = xformrot( X);
		}	
		endeff(i) = rel2abs( M, M.getN()-1, Q(i), RR );	
	}
	
	/* initialize IL */
	ilInit();
	
    /* wait for events */ 
    while ( !done )
	{
	    /* handle the events in the queue */

	    while ( SDL_PollEvent( &event ) )
		{
		    switch( event.type )
			{
			case SDL_ACTIVEEVENT:
			    /* Something's happend with our focus
			     * If we lost focus or we are iconified, we
			     * shouldn't draw the screen
			     */
			    if ( event.active.gain == 0 )
				{
					isActive = FALSE;
					timeAtInactive = timer.toc();
				}
			    else
				{
					isActive = TRUE;
					durationInactive += timer.toc() - timeAtInactive;
				}
			    break;			    
			case SDL_VIDEORESIZE:
			    /* handle resize event */
			    surface = SDL_SetVideoMode( event.resize.w,
							event.resize.h,
							16, videoFlags );
			    if ( !surface )
				{
				    fprintf( stderr, "Could not get a surface after resize: %s\n", SDL_GetError( ) );
				    Quit( 1 );
				}
			    resizeWindow( event.resize.w, event.resize.h );
			    break;
			case SDL_KEYDOWN:
			    /* handle key presses */
			    handleKeyPress( &event.key.keysym );
			    break;
			case SDL_QUIT:
			    /* handle quit requests */
			    done = TRUE;
			    break;
			default:
			    break;
			}
		}

		durationActive = timer.toc() - durationInactive;

	    /* draw the scene */
	    if ( isActive )
		{
			/* Clear The Screen And The Depth Buffer */
			glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

			drawGLScene( M, pos, axis, endeff); // YAYYYYYYYY
			if ( saveImgs && tidx != maxTidx )
			{
				sprintf(imgname,"anim%d.tga",tidx);
				Screendump( imgname, SCREEN_WIDTH, SCREEN_HEIGHT);
			   			    
			   	ilLoadImage(imgname);
			   	sprintf(imgname,"anim%d.png",tidx);
			   	ilSave(IL_PNG,imgname);
		   	}
		}
	}

    /* clean ourselves up and exit */
    Quit( 0 );

    /* Should never get here */
    return( 0 );
}

