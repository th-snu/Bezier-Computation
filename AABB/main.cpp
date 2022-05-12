#include <GL/glut.h>
#include <stdio.h>
#include <math.h>
#include "biarc_approx.h"
#include "aabb.h"

#define RES 100

CubicBezierCurve curve;
GLsizei width = 1280, height = 960;
int edit_ctrlpts_idx = -1;
bool isDrawControlMesh = true;
bool isDottedLine = false;
bool isDrawBiarcs = false;
bool isDrawAABB = true;
bool isDrawHierarchy = false;
int subdivision_power = 4;

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

REAL get_AABB(const CubicBezierCurve &seg, const Arc arc, AABB& box, bool isDraw){
	REAL error;

	glColor3ub(64, 192, 0);
	if (arc.radius == arc.radius && (arc.end - arc.begin > M_PI / 2.0 || arc.end - arc.begin <= 0.001)){
		// replace wide-angled arc with a line instead
		// extremely narrow-angled arc has extreme error, so use line instead of them as well since they're close to a line

		const CubicBezierCurve &dseg = seg;
		CubicBezierCurve approx;
		copy_point(dseg.control_pts[0], approx.control_pts[0]);
		copy_point(dseg.control_pts[3], approx.control_pts[3]);
		division_point(dseg.control_pts[0], dseg.control_pts[3], 1.0/3.0, approx.control_pts[1]);
		division_point(dseg.control_pts[0], dseg.control_pts[3], 2.0/3.0, approx.control_pts[2]);
		error = bezier_error_bound(&dseg, &approx);

		box.x[0] = std::min(approx.control_pts[0][0], approx.control_pts[3][0]);
		box.x[1] = std::max(approx.control_pts[0][0], approx.control_pts[3][0]);
		box.y[0] = std::min(approx.control_pts[0][1], approx.control_pts[3][1]);
		box.y[1] = std::max(approx.control_pts[0][1], approx.control_pts[3][1]);

		if (isDraw){
			glBegin(GL_LINE_STRIP);
			glVertex2f(approx.control_pts[0][0], approx.control_pts[0][1]);
			glVertex2f(approx.control_pts[3][0], approx.control_pts[3][1]);
			glEnd();
		}
	}
	else {
		box = get_arc_aabb(&arc);

		error = arc_approx_error_bound(&arc, &seg);

		if (isDraw){
			if (arc.radius != arc.radius){
				glBegin(GL_LINE_STRIP);
				glVertex2f(arc.center[0], arc.center[1]);
				glVertex2f(arc.begin, arc.end);
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

	return error;
}

void draw_AABB(const AABB &box){
	glBegin(GL_LINE_STRIP);
	glVertex2f(box.x[0], box.y[0]);
	glVertex2f(box.x[0], box.y[1]);
	glVertex2f(box.x[1], box.y[1]);
	glVertex2f(box.x[1], box.y[0]);
	glVertex2f(box.x[0], box.y[0]);
	glEnd();
}

void build_hierarchy(std::shared_ptr<Hierarchy> h, int power){
	std::vector<CubicBezierCurve> curves;
	CubicBezierCurve &seg = h->curve;

	Point inflect;
	Arc arc1, arc2;
	get_biarc_inflect(&seg, inflect);
	// draw_arc_center(&seg, inflect, &arc1, &arc2);
	to_biarc(&seg, inflect, &arc1, &arc2);

	// divide seg at inflection point
	std::vector<CubicBezierCurve> segs;
	subdivide(&seg, segs, 1);

	auto leftH = std::make_shared<Hierarchy>();
	auto rightH = std::make_shared<Hierarchy>();
	leftH->curve = segs[0];
	rightH->curve = segs[1];

	if (power == 0){
		REAL error1 = get_AABB(segs[0], arc1, leftH->box, isDrawBiarcs && !power);
		REAL error2 = get_AABB(segs[1], arc2, rightH->box, isDrawBiarcs && !power);

		leftH->box.x[0] -= error1;
		leftH->box.x[1] += error1;
		leftH->box.y[0] -= error1;
		leftH->box.y[1] += error1;
		rightH->box.x[0] -= error2;
		rightH->box.x[1] += error2;
		rightH->box.y[0] -= error2;
		rightH->box.y[1] += error2;
	}
	else {
		build_hierarchy(leftH, power - 1);
		build_hierarchy(rightH, power - 1);
	}
	h->box = combine(leftH->box, rightH->box);
	h->left = leftH;
	h->right = rightH;

	if (isDrawAABB && (!power || isDrawHierarchy)){
		if (power % 2 == 0)
			glColor3ub(64, 0, 255);
		else glColor3ub(255, 0, 64);
		draw_AABB(h->box);
	}
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

	auto root = std::make_shared<Hierarchy>();

	root->curve = curve;
	build_hierarchy(root, subdivision_power);

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
	case '1':
		isDrawBiarcs ^= true;
		break;
	case '2':
		isDrawAABB ^= true;
		break;
	case '3':
		isDrawHierarchy ^= true;
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
	glutCreateWindow("AABB Hierarchy");

	init();
	glutReshapeFunc(reshape_callback);
	glutMouseFunc(mouse_callback);
	glutMotionFunc(mouse_move_callback);
	glutDisplayFunc(display_callback);
	glutKeyboardFunc(keyboard_callback);
	glutMainLoop();

	return 0;
}
