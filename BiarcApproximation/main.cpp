#include <GL/glut.h>
#include <stdio.h>
#include <math.h>
#include "curve.h"

#define RES 100

CubicBezierCurve curve;
GLsizei width = 1280, height = 960;
int edit_ctrlpts_idx = -1;
bool isDrawControlMesh = true;
bool isDottedLine = false;
bool isDrawSegControl = false;
bool isDrawInflectionPoint = false;
bool isDrawTangentLine = false;
bool isDrawBisectionCircle = false;
bool isDrawBiarcs = true;
int subdivision_power = 2;
bool usePreciseInflection = false;

int hit_index(CubicBezierCurve *curve, int x, int y)
{
	for (int i = 0; i < 4; i++)
	{
		REAL tx = curve->control_pts[i][0] - x;
		REAL ty = curve->control_pts[i][1] - y;
		if ((tx * tx + ty * ty) < 30) return i;
	}
	return -1;
}

void init()
{
	SET_PT2(curve.control_pts[0], 200, 300);
	SET_PT2(curve.control_pts[1], 450, 700);
	SET_PT2(curve.control_pts[2], 750, 700);
	SET_PT2(curve.control_pts[3], 1000, 300);

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, width, 0, height);
}

void reshape_callback(GLint nw, GLint nh)
{
	width = nw;
	height = nh;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, width, 0, height);
}

void draw_bezier_seg(std::vector<CubicBezierCurve> &curves){
	glColor3ub(0, 0, 0);
	for (CubicBezierCurve seg : curves){
		glBegin(GL_LINES);
		for (int i = 0; i <= RES; i++)
		{
			Point pt;
			const REAL t = (REAL)i / (REAL)RES;
			evaluate(&curve, t, pt);
			glVertex2f(pt[0], pt[1]);
		}
		glEnd();
	}
}

void draw_line(const Point p, const Point tan){
	Point pt;
	glBegin(GL_LINE_STRIP);
	pt[0] = p[0] - tan[0] * 10000.0;
	pt[1] = p[1] - tan[1] * 10000.0;
	glVertex2f(pt[0], pt[1]);
	pt[0] = p[0] + tan[0] * 10000.0;
	pt[1] = p[1] + tan[1] * 10000.0;
	glVertex2f(pt[0], pt[1]);
	glEnd();
}

void draw_tangent_line(const std::vector<CubicBezierCurve> &curves){
	glColor3ub(0, 255, 0);
	for (auto seg: curves){
		Point tan_begin, tan_end;
		get_tangent(&seg, tan_begin, tan_end);

		Point pt;
		draw_line(seg.control_pts[0], tan_begin);
		draw_line(seg.control_pts[3], tan_end);
	}
}

void draw_segment_controls(const std::vector<CubicBezierCurve> &curves){
	for (auto seg: curves){
		/* segment control mesh */
		if (isDrawControlMesh)
		{
			glColor3ub(128, 0, 0);
			glBegin(GL_LINE_STRIP);
			for (int i = 0; i < 4; i++)
			{
				REAL *pt = seg.control_pts[i];
				glVertex2f(pt[0], pt[1]);
			}
			glEnd();
		}

		/* segment control points */
		glColor3ub(0, 0, 128);
		glPointSize(5.0);
		glBegin(GL_POINTS);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = seg.control_pts[i];
			glVertex2f(pt[0], pt[1]);
		}
		glEnd();
		glutSwapBuffers();
	}
}

void draw_bisection_circle(const std::vector<CubicBezierCurve> &curves){
	glColor3ub(255, 255, 0);
	for (CubicBezierCurve seg: curves){
		Point l1_tan, l1_center, l2_tan, l2_center, center;
		get_line_segs(&seg, l1_tan, l1_center, l2_tan, l2_center);

		// draw_line(l1_center, l1_tan);
		// draw_line(l2_center, l2_tan);
		get_line_intersection(l1_tan, l1_center, l2_tan, l2_center, center);
		glPointSize(10.0);
		glBegin(GL_POINTS);
		// glVertex2f(l1_center[0], l1_center[1]);
		// glVertex2f(l2_center[0], l2_center[1]);
		glVertex2f(center[0], center[1]);
		glEnd();

		glBegin(GL_LINE_STRIP);
		for (int i = 0; i <= RES; i++)
		{
			Point pt;
			Point r_vec;
			subtract_point(center, seg.control_pts[0], r_vec);
			REAL radius = norm(r_vec);

			const REAL angle = (REAL)M_PI * 2 * i / (REAL)RES;
			pt[0] = center[0] + radius * cos(angle);
			pt[1] = center[1] + radius * sin(angle);
			glVertex2f(pt[0], pt[1]);
		}
		glEnd();
	}
}

void draw_inflect(const std::vector<CubicBezierCurve> &curves)
{
	glColor3ub(255, 0, 255);
	glPointSize(5.0);
	for (CubicBezierCurve seg: curves){
		Point inflect;
		Point l1_tan, l1_center, l2_tan, l2_center;
		get_line_segs(&seg, l1_tan, l1_center, l2_tan, l2_center);
		get_biarc_inflect(&seg, inflect, usePreciseInflection);
		glBegin(GL_POINTS);
		glVertex2f(inflect[0], inflect[1]);
		glEnd();
	}
}

void draw_arc_center(const CubicBezierCurve *curve, const Point &inflect, Arc *arc1, Arc *arc2)
{
	const Point &arc1_point = curve->control_pts[0];
	const Point &arc2_point = curve->control_pts[3];
	set_arc_center(curve, arc1_point, arc2_point, inflect, arc1, arc2);

	glColor3ub(0, 255, 255);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	glVertex2f(arc1->center[0], arc1->center[1]);
	glVertex2f(arc2->center[0], arc2->center[1]);
	glEnd();
}

void display_callback()
{
	/* curve */
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3ub(0, 0, 0);
	if (isDottedLine)
		glBegin(GL_LINES);
	else
		glBegin(GL_LINE_STRIP);	
	for (int i = 0; i <= RES; i++)
	{
		Point pt;
		const REAL t = (REAL)i / (REAL)RES;
		evaluate(&curve, t, pt);
		glVertex2f(pt[0], pt[1]);
	}
	glEnd();

	/* subdivision */
	std::vector<CubicBezierCurve> curves;

	subdivide(&curve, curves, subdivision_power);
	// draw_bezier_seg(curves);
	if (isDrawSegControl) draw_segment_controls(curves);
	if (isDrawTangentLine) draw_tangent_line(curves);
	if (isDrawBisectionCircle) draw_bisection_circle(curves);
	if (isDrawInflectionPoint) draw_inflect(curves);

	/* biarcs */
	if (isDrawBiarcs){
		glColor3ub(0, 255, 255);
		std::vector<Arc> arcs;
		for (CubicBezierCurve seg : curves){
			Point inflect;
			Arc arc1, arc2;
			get_biarc_inflect(&seg, inflect, usePreciseInflection);
			// draw_arc_center(&seg, inflect, &arc1, &arc2);
			to_biarc(&seg, inflect, &arc1, &arc2);
			arcs.push_back(arc1);
			arcs.push_back(arc2);
		}
		for (auto arc: arcs){
			if (arc.radius != arc.radius){
				glBegin(GL_LINE_STRIP);
				glVertex2f(arc.center[0], arc.center[1]);
				glVertex2f(arc.center[0] + arc.begin, arc.center[1] + arc.end);
				glEnd();
			}
			else {
				glBegin(GL_LINE_STRIP);
				for (int i = 0; i <= RES; i++)
				{
					Point pt;
					const REAL t = (REAL)i / (REAL)RES;
					const REAL angle = t * (arc.end - arc.begin) + arc.begin;
					pt[0] = arc.center[0] + arc.radius * cos(angle);
					pt[1] = arc.center[1] + arc.radius * sin(angle);
					glVertex2f(pt[0], pt[1]);
				}
				glEnd();
			}
		}
	}

	/* control mesh */
	if (isDrawControlMesh)
	{
		glColor3ub(255, 0, 0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = curve.control_pts[i];
			glVertex2f(pt[0], pt[1]);
		}
		glEnd();
	}

	/* control pts */
	glColor3ub(0, 0, 255);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < 4; i++)
	{
		REAL *pt = curve.control_pts[i];
		glVertex2f(pt[0], pt[1]);
	}
	glEnd();
	glutSwapBuffers();
}

// void glutMouseFunc(void (*func)(int button, int state, int x, int y));
void mouse_callback(GLint button, GLint action, GLint x, GLint y)
{
	if (GLUT_LEFT_BUTTON == button)
	{
		switch (action)
		{
			case GLUT_DOWN:
				edit_ctrlpts_idx = hit_index(&curve, x, height - y);
				break;
			case GLUT_UP:
				edit_ctrlpts_idx = -1;
				break;
			default: break;
		}
	}
	glutPostRedisplay();
}

// void glutMotionFunc(void (*func)(int x, int y));
void mouse_move_callback(GLint x, GLint y)
{
	if (edit_ctrlpts_idx != -1)
	{
		curve.control_pts[edit_ctrlpts_idx][0] = x;
		curve.control_pts[edit_ctrlpts_idx][1] = height - y;
	}
	glutPostRedisplay();
}

// void glutKeyboardFunc(void (*func)(unsigned char key, int x, int y));
void keyboard_callback(unsigned char key, int x, int y)
{
	switch(key)
	{
	case 'i': case 'I':
		SET_PT2(curve.control_pts[0], 200, 300);
		SET_PT2(curve.control_pts[1], 450, 700);
		SET_PT2(curve.control_pts[2], 750, 700);
		SET_PT2(curve.control_pts[3], 1000, 300);
		break;
	case 'l': case 'L':
		isDottedLine ^= true;
		break;
	case 'c': case 'C':
		isDrawControlMesh ^= true;
		break;
	// case 'a': case 'A':
	// 	usePreciseInflection ^= true;
	// 	break;
	case '1':
		isDrawSegControl ^= true;
		break;
	case '2':
		isDrawTangentLine ^= true;
		break;
	case '3':
		isDrawBisectionCircle ^= true;
		break;
	case '4':
		isDrawInflectionPoint ^= true;
		break;
	case '5':
		isDrawBiarcs ^= true;
		break;
	case '=': case '+':
		subdivision_power += 1;
		break;
	case '-':
		subdivision_power -= 1;
		if (subdivision_power < 0) subdivision_power = 0;
		break;
	case (27): exit(0); break;
	default: break;
	}
	glutPostRedisplay();
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(width, height);
	glutCreateWindow("Bezier Editor");

	init();
	glutReshapeFunc(reshape_callback);
	glutMouseFunc(mouse_callback);
	glutMotionFunc(mouse_move_callback);
	glutDisplayFunc(display_callback);
	glutKeyboardFunc(keyboard_callback);
	glutMainLoop();

	return 0;
}
