#include "visual3d.h"

using namespace std;
using namespace Eigen;

MatGraph3d* matGraph3d = NULL;
void SetMatGraph3d(MatGraph3d *_matGraph3d){
	matGraph3d = _matGraph3d;
}

GLfloat camerapos[3] = {-4, 6, -4};
GLfloat cameradir[3] = {sqrt(2)/2, -0.4, sqrt(2)/2};
const GLfloat Pi = 3.1415926;


void InitGL(int argc, char** argv){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(500, 500);
	glutCreateWindow("OpenGL 3D View");
	
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.3, 0.3, 0.3, 0.0);
	//glMatrixMode(GL_PROJECTION);

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutMainLoop();
	return;
}

void display()
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(camerapos[0],camerapos[1],camerapos[2], 
		camerapos[0]+cameradir[0],camerapos[1]+cameradir[1],camerapos[2] + cameradir[2],
		0,1,0);

	//glEnable(GL_DEPTH_TEST);
	
	

	if(matGraph3d)
		matGraph3d->ShowResult();
	glutSwapBuffers();
}


float movestep = 0.1f;
void keyboard(unsigned char key, int x, int y)
{
	switch(key){
	case 'w':
		camerapos[0]+=movestep*cameradir[0];
		camerapos[1]+=movestep*cameradir[1];
		camerapos[2]+=movestep*cameradir[2];
		break;
	case 's':
		camerapos[0]-=movestep*cameradir[0];
		camerapos[1]-=movestep*cameradir[1];
		camerapos[2]-=movestep*cameradir[2];
		break;
	case 'd':
		camerapos[0]-=movestep*cameradir[2];
		//camerapos[1]+=movestep*cameradir[1];
		camerapos[2]+=movestep*cameradir[0];
		break;
	case 'a':
		camerapos[0]+=movestep*cameradir[2];
		//camerapos[1]+=movestep*cameradir[1];
		camerapos[2]-=movestep*cameradir[0];
		break;
	case 'e':
		camerapos[1]+=movestep;
		break;
	case 'q':
		camerapos[1]-=movestep;
		break;
	}
	glutPostRedisplay();
}

int mx,my = 0;
bool click = false;
void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_DOWN) {
		mx = x;
		my = y;
		click = true; 
	}
	else
		click = false;
}

GLfloat rotArr[2] = {Pi/4, 0};//2pi pi

void _updateCamDir(){
	cameradir[1] = sin(rotArr[1]);
	cameradir[0] = cos(rotArr[1])*sin(rotArr[0]);
	cameradir[2] = cos(rotArr[1])*cos(rotArr[0]);
}

void motion(int x, int y)
{
	if (click)	
	{

		rotArr[1] += (y - my) / 40.0f;
		rotArr[1] = rotArr[1]>Pi/2-0.1 ? Pi/2:rotArr[1];
		rotArr[1] = rotArr[1]<-Pi/2+0.1 ? -Pi/2:rotArr[1];
		
		rotArr[0] += (x - mx) / 40.0f;
		mx = x;
		my = y;
	}
	_updateCamDir();
	glutPostRedisplay();
}

void reshape(GLsizei w, GLsizei h)
{
	if (h == 0)
		h = 1;
	float aspect = (float)w / h;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(75, aspect, 0.1, 100);

}


double InterpolateZero2One(double min, double max, double input){
	if(min == max)
		return input;
	return (input - min)/(max - min);
}

void SetVertexColor(double weight){
	if(weight > 1 || weight < 0)
	{
		glColor3f(1.0f,1.0f,1.0f);
		return;
	}
	else{
		glColor4f(weight, 0.0f, 1.0f - weight, .9f);		
	}
}

void SetMatVertexColor(double minHeight, double maxHeight, double height){
	SetVertexColor(InterpolateZero2One(minHeight, maxHeight, height));
}


void DrawCube(double x, double y, double z, double len){
	glBegin(GL_QUADS);    
    glNormal3f( 0.0F, 0.0F, 1.0F);    
    glVertex3f(x+0.5f*len,y+0.5f*len,z+0.5f*len);    
    glVertex3f(x-0.5f*len,y+0.5f*len,z+0.5f*len);    
    glVertex3f(x-0.5f*len,y-0.5f*len,z+0.5f*len);    
    glVertex3f(x+0.5f*len,y-0.5f*len,z+0.5f*len);    
    //1----------------------------    
    glNormal3f( 0.0F, 0.0F,-1.0F);    
    glVertex3f(x-0.5f*len,y-0.5f*len,z-0.5f*len);    
    glVertex3f(x-0.5f*len,y+0.5f*len,z-0.5f*len);    
    glVertex3f(x+0.5f*len,y+0.5f*len,z-0.5f*len);    
    glVertex3f(x+0.5f*len,y-0.5f*len,z-0.5f*len);    
    //2----------------------------    
    glNormal3f( 0.0F, 1.0F, 0.0F);    
    glVertex3f(x+0.5f*len,y+0.5f*len,z+0.5f*len);    
    glVertex3f(x+0.5f*len,y+0.5f*len,z-0.5f*len);    
    glVertex3f(x-0.5f*len,y+0.5f*len,z-0.5f*len);    
    glVertex3f(x-0.5f*len,y+0.5f*len,z+0.5f*len);    
    //3----------------------------    
    glNormal3f( 0.0F,-1.0F, 0.0F);    
    glVertex3f(x-0.5f*len,y-0.5f*len,z-0.5f*len);    
    glVertex3f(x+0.5f*len,y-0.5f*len,z-0.5f*len);    
    glVertex3f(x+0.5f*len,y-0.5f*len,z+0.5f*len);    
    glVertex3f(x-0.5f*len,y-0.5f*len,z+0.5f*len);    
    //4----------------------------    
    glNormal3f( 1.0F, 0.0F, 0.0F);    
    glVertex3f(x+0.5f*len,y+0.5f*len,z+0.5f*len);    
    glVertex3f(x+0.5f*len,y-0.5f*len,z+0.5f*len);    
    glVertex3f(x+0.5f*len,y-0.5f*len,z-0.5f*len);    
    glVertex3f(x+0.5f*len,y+0.5f*len,z-0.5f*len);    
    //5----------------------------    
    glNormal3f(-1.0F, 0.0F, 0.0F);
    glVertex3f(x-0.5f*len,y-0.5f*len,z-0.5f*len);    
    glVertex3f(x-0.5f*len,y-0.5f*len,z+0.5f*len);    
    glVertex3f(x-0.5f*len,y+0.5f*len,z+0.5f*len);    
    glVertex3f(x-0.5f*len,y+0.5f*len,z-0.5f*len);    
    //6----------------------------*/    
    glEnd();
}

void MatGraph3d::SetMask(int *mask){
	_mask = mask;
}


void MatGraph3d::ShowResult(){
	int X = 0, Y = 0, Z=0;//设置循环变量
	double x, y, z;

	//draw grid
	glColor4f(1.0f, 1.0f, 1.0f, 0.8f);// 重置颜色
	glBegin( GL_LINES );
	for ( X = 0; X < _dimx-1; X ++ )
	for ( Z = 0; Z < _dimz-1; Z ++ )
	{
		glVertex3f(_scale * X    , 0, _scale * Z    );
		glVertex3f(_scale * X    , 0, _scale * (Z+1));

		glVertex3f(_scale * X    , 0, _scale * (Z+1));
		glVertex3f(_scale * (X+1), 0, _scale * (Z+1));

		glVertex3f(_scale * (X+1), 0, _scale * (Z+1));
		glVertex3f(_scale * (X+1), 0, _scale * Z    );

		glVertex3f(_scale * (X+1), 0, _scale * Z    );
		glVertex3f(_scale * X    , 0, _scale * Z    );
	}
	glColor3f(1,0,0);
	glVertex3f(_scale * 2 * _dimx, 0, 0);
	glVertex3f(-_scale * 2 * _dimx, 0, 0);
	
	glColor3f(0,1,0);
	glVertex3f(0, 0, _scale * 2 * _dimz);
	glVertex3f(0, 0, -_scale * 2 * _dimz);

	glEnd();

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	//glBlendEquation(GL_FUNC_ADD);


	double maxHeight = -1000000, minHeight = 1000000;
	for(int i = 0; i<_dimx; i++)
	for(int j = 0; j<_dimy; j++)
	for(int k = 0; k<_dimz; k++)
	{
		maxHeight = maxHeight>_mat[k*_dimx*_dimy+j*_dimx+i] ? maxHeight:_mat[k*_dimx*_dimy+j*_dimx+i];
		minHeight = minHeight<_mat[k*_dimx*_dimy+j*_dimx+i] ? minHeight:_mat[k*_dimx*_dimy+j*_dimx+i];
	}
	// if(_renderType)//选择渲染模式
	// 	glBegin( GL_QUADS );// 渲染为四边形
	// else
	// 	glBegin( GL_LINES );// 渲染为直线

	for ( X = 0; X < _dimx; X ++ )
	for ( Y = 0; Y < _dimy; Y ++ )
	for ( Z = 0; Z < _dimz; Z ++ )
	{
		
		// 绘制(x,y)处的顶点
		// 获得(x,y,z)坐标
		x = _scale * X;
		y = _scale * Y;
		z = _scale * Z;
		if(_mask!=NULL && (_mask[Z*_dimx*_dimy+Y*_dimx+X]!=0))
			glColor3f(1.0f,1.0f,1.0f);
		else
			SetMatVertexColor(minHeight, maxHeight, _mat[Z*_dimx*_dimy+Y*_dimx+X]);

		DrawCube(x,y,z,_scale/6);
	}
	//glEnd();


	glDisable(GL_BLEND);
}
