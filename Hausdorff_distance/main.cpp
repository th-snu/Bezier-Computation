#include <GL/glut.h>
#include <stdio.h>
#include <math.h>
#include <random>
#include "hausdorff.h"

#define RES 100
#define NUM_SAMPLES 10

CubicBezierCurve curve1;
CubicBezierCurve curve2;
GLsizei width = 1280, height = 960;
float offset_x = 0, offset_y = 0;
double mag = 1.0;
int edit_ctrlpts_idx = -1;
int edit_curve_idx = -1;
bool isDrawControlMesh = true;
bool isDottedLine = false;
bool isDrawBiarcs = false;
int subdivision_power = 6;
int text_line = 0;
int old_x, old_y;

int hit_index(CubicBezierCurve *curve, int x, int y)
{
	float fx = x * mag;
	float fy = y * mag;
	fx += offset_x;
	fy += offset_y;
	for (int i = 0; i < 4; i++)
	{
		REAL tx = curve->control_pts[i][0] - fx;
		REAL ty = curve->control_pts[i][1] - fy;
		if ((tx * tx + ty * ty) < 30.0 * mag * mag) return i;
	}
	return -1;
}

void init_points()
{
	SET_PT2(curve1.control_pts[0], 200, 300);
	SET_PT2(curve1.control_pts[1], 450, 700);
	SET_PT2(curve1.control_pts[2], 750, 700);
	SET_PT2(curve1.control_pts[3], 1000, 300);

	SET_PT2(curve2.control_pts[0], 200, 700);
	SET_PT2(curve2.control_pts[1], 450, 300);
	SET_PT2(curve2.control_pts[2], 750, 300);
	SET_PT2(curve2.control_pts[3], 1000, 700);
	mag = 1.0;
	offset_x = 0;
	offset_y = 0;
}

void init()
{
	init_points();
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
	gluOrtho2D(offset_x, offset_x + width * mag, offset_y, offset_y + height * mag);
}

void draw_curve(CubicBezierCurve& curve){
	if (isDottedLine)
		glBegin(GL_LINES);
	else
		glBegin(GL_LINE_STRIP);	
	for (int i = 0; i <= RES; i++){
		Point pt;
		const REAL t = (REAL)i / (REAL)RES;
		evaluate(&curve, t, pt);
		glVertex2f(pt[0], pt[1]);
	}
	glEnd();
}

void draw_text(std::string str){
	glColor3ub(0, 0, 0);
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D( 0, width, 0, height );

	glMatrixMode( GL_MODELVIEW );
	glPushMatrix();
	glLoadIdentity();
	glRasterPos2i( 32, height - 32 - 32 * text_line );  // move in 10 pixels from the left and bottom edges
	for ( int i = 0; i < str.size(); ++i ) {
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, str[i]);
	}
	glPopMatrix();

	glMatrixMode( GL_PROJECTION );
	glPopMatrix();
	glMatrixMode( GL_MODELVIEW );
	text_line += 1;
}

typedef std::pair<REAL, std::pair<std::pair<REAL, REAL>, bool>> max_pair;
void draw_hausdorff_distance(const CubicBezierCurve &curve1, const CubicBezierCurve &curve2){
	std::priority_queue<max_pair> q;
	
	Point tmp_pt1, tmp_pt2, tmp_pt3, tmp_pt4;
	REAL t1_lower_bound = sample_lower_bound(curve1, curve2, NUM_SAMPLES, tmp_pt1, tmp_pt2);
	REAL t2_lower_bound = sample_lower_bound(curve2, curve1, NUM_SAMPLES, tmp_pt3, tmp_pt4);
	REAL lower_bound = std::max(t1_lower_bound, t2_lower_bound);
	REAL upper_bound = bezier_error_bound(&curve1, &curve2);
	Point bound1, bound2;
	if (t1_lower_bound > t2_lower_bound){
		copy_point(tmp_pt1, bound1);
		copy_point(tmp_pt2, bound2);
	}
	else{
		copy_point(tmp_pt3, bound1);
		copy_point(tmp_pt4, bound2);
	}

	// Bezier segment from tree1 will be marked with False, and bezier segment from tree2 will be marked with True
	q.push(std::make_pair(upper_bound, std::make_pair(std::make_pair(0.0, 1.0), false)));
	q.push(std::make_pair(upper_bound, std::make_pair(std::make_pair(0.0, 1.0), true)));

	int last_update = 0;

	while (!q.empty()){
		REAL curr_bound = q.top().first;
		if (curr_bound < lower_bound + 1e-5){
			upper_bound = curr_bound;
			break;
		}
		if (upper_bound > curr_bound) last_update = 0;
		else last_update += 1;
		if (last_update > 100) break;
		upper_bound = curr_bound;

		auto local_t = q.top().second.first;
		auto idx = q.top().second.second;
		auto proj_target = idx ? curve1 : curve2;
		auto local_curve = idx ? curve2 : curve1;
		q.pop();

		CubicBezierCurve local_seg = subcurve_by_endpoint(local_curve, local_t.first, local_t.second);

		CubicBezierCurve left_seg, right_seg;
		subdivide(&local_seg, &left_seg, &right_seg);

		// Add child nodes to priority queue, compute upperbound by projecting two end points to other bezier
		REAL t1 = projection(left_seg.control_pts[0], proj_target);
		REAL t2 = projection(left_seg.control_pts[3], proj_target);
		REAL t3 = projection(right_seg.control_pts[3], proj_target);

		CubicBezierCurve c1 = subcurve_by_endpoint(proj_target, std::min(t1, t2), std::max(t1, t2));
		CubicBezierCurve c2 = subcurve_by_endpoint(proj_target, std::min(t2, t3), std::max(t3, t2));

		REAL upper_bound_left = bezier_error_bound(&(left_seg), &c1);
		upper_bound_left = std::min(upper_bound_left, curr_bound);
		REAL upper_bound_right = bezier_error_bound(&(right_seg), &c2);
		upper_bound_right = std::min(upper_bound_right, curr_bound);

		Point sample1, sample2;
		REAL lower_bound_left = sample_lower_bound(left_seg, proj_target, NUM_SAMPLES, sample1, sample2);
		if (lower_bound < lower_bound_left) {
			lower_bound = lower_bound_left;
			copy_point(sample1, bound1);
			copy_point(sample2, bound2);
		}
		REAL lower_bound_right = sample_lower_bound(right_seg, proj_target, NUM_SAMPLES, sample1, sample2);
		if (lower_bound < lower_bound_right) {
			lower_bound = lower_bound_right;
			copy_point(sample1, bound1);
			copy_point(sample2, bound2);
		}

		q.push(std::make_pair(upper_bound_left, std::make_pair(std::make_pair(local_t.first, (local_t.first + local_t.second) / 2.0), idx)));
		q.push(std::make_pair(upper_bound_right, std::make_pair(std::make_pair((local_t.first + local_t.second) / 2.0, local_t.second), idx)));
	}
	
	std::string distance = "Distance: " + std::to_string((upper_bound + lower_bound) / 2);
	std::string error = "Error: " + std::to_string((upper_bound - lower_bound) / 2.0);
	glColor3ub(255, 0, 0);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	glVertex2f(bound1[0], bound1[1]);
	glPointSize(1.0);
	glEnd();
	glLineWidth(5.0);
	glBegin(GL_LINES);
	glVertex2f(bound1[0], bound1[1]);
	glVertex2f(bound2[0], bound2[1]);
	glEnd();
	glLineWidth(1.0);
	draw_text(distance);
	draw_text(error);
}

void display_callback()
{
	text_line = 0;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(offset_x, offset_x + width * mag, offset_y, offset_y + height * mag);

	/* curve */
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3ub(0, 0, 0);

	glLineWidth(2.0);
	draw_curve(curve1);
	draw_curve(curve2);
	glLineWidth(1.0);

	/* control pts */
	glColor3ub(0, 0, 255);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < 4; i++)
	{
		REAL *pt = curve1.control_pts[i];
		glVertex2f(pt[0], pt[1]);
	}
	for (int i = 0; i < 4; i++)
	{
		REAL *pt = curve2.control_pts[i];
		glVertex2f(pt[0], pt[1]);
	}
	glEnd();

	REAL t = projection(curve1.control_pts[0], curve2);
	Point proj;
	evaluate(&curve2, t, proj);
	glBegin(GL_LINES);
	glVertex2f(curve1.control_pts[0][0], curve1.control_pts[0][1]);
	glVertex2f(proj[0], proj[1]);
	glEnd();
	// draw_hausdorff_distance(curve1, curve2);

	draw_text("Magnification: " + std::to_string(1.0/mag));

	/* control mesh */
	if (isDrawControlMesh)
	{
		glColor3ub(255, 0, 0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = curve1.control_pts[i];
			glVertex2f(pt[0], pt[1]);
		}
		glEnd();
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = curve2.control_pts[i];
			glVertex2f(pt[0], pt[1]);
		}
		glEnd();
	}
	glutSwapBuffers();
}

void update_mag(double val){
	double old_mag = mag;
	double old_offset_x = offset_x;
	double old_offset_y = offset_y;
	mag *= val;
	offset_x = old_offset_x + (width * old_mag - width * mag) / 2;
	offset_y = old_offset_y + (height * old_mag - height * mag) / 2;
}

// void glutMouseFunc(void (*func)(int button, int state, int x, int y));
void mouse_callback(GLint button, GLint action, GLint x, GLint y)
{
	if (button == 3 || button == 4){
		if (action == GLUT_UP) return;
		else {
			if (button == 3) update_mag(1.0/1.1);
			else update_mag(1.1);
		}
	}
	if (GLUT_LEFT_BUTTON == button)
	{
		switch (action)
		{
			case GLUT_DOWN:
				{
					edit_ctrlpts_idx = hit_index(&curve1, x, height - y);
					edit_curve_idx = 1;
				}
				if (edit_ctrlpts_idx == -1){
					edit_ctrlpts_idx = hit_index(&curve2, x, height - y);
					edit_curve_idx = 2;
				}
				if (edit_ctrlpts_idx == -1){
					edit_curve_idx = -1;
				}
				break;
			case GLUT_UP:
				edit_ctrlpts_idx = -1;
				edit_curve_idx = -1;
				break;
			default: break;
		}
	}
	glutPostRedisplay();
	old_x = x;
	old_y = y;
}

// void glutMotionFunc(void (*func)(int x, int y));
void mouse_move_callback(GLint x, GLint y)
{
	float fy = height - y;
	float fx = x * mag;
	fy *= mag;
	fx += offset_x;
	fy += offset_y;
	if (edit_ctrlpts_idx != -1)
	{
		auto& curve = edit_curve_idx == 1 ? curve1 : curve2;
		curve.control_pts[edit_ctrlpts_idx][0] = fx;
		curve.control_pts[edit_ctrlpts_idx][1] = fy;
	}
	else {
		offset_x -= (float)(x - old_x) * mag;
		offset_y -= (float)(old_y - y) * mag;
	}
	glutPostRedisplay();
	old_x = x;
	old_y = y;
}

void snapshot_points(){
	std::ofstream fout("points.txt");
	for (int i = 0; i < 4; i++){
		fout << curve1.control_pts[i][0] << " " << curve1.control_pts[i][1] << std::endl;
	}
	for (int i = 0; i < 4; i++){
		fout << curve2.control_pts[i][0] << " " << curve2.control_pts[i][1] << std::endl;
	}
	fout.close();
}

void load_points(){
	try{
		std::ifstream fin("points.txt");
		for (int i = 0; i < 4; i++){
			fin >> curve1.control_pts[i][0] >> curve1.control_pts[i][1];
		}
		for (int i = 0; i < 4; i++){
			fin >> curve2.control_pts[i][0] >> curve2.control_pts[i][1];
		}
		fin.close();
	}
	catch (...){
		std::cout << "No points file found.\n";
	}
}

// void glutKeyboardFunc(void (*func)(unsigned char key, int x, int y));
void keyboard_callback(unsigned char key, int x, int y)
{
	switch(key)
	{
	case 'i': case 'I':
		init_points();
		break;
	case 'l': case 'L':
		isDottedLine ^= true;
		break;
	case 'c': case 'C':
		isDrawControlMesh ^= true;
		break;
	case '1':
		snapshot_points();
		break;
	case '2':
	    load_points();
		break;
	case '=': case '+':
		subdivision_power += 1;
		break;
	case '-':
		subdivision_power -= 1;
		if (subdivision_power < 0) subdivision_power = 0;
		break;
	case ',':
		update_mag(1.5);
		break;
	case '.':
		update_mag(1.0/1.5);
		break;
	case 'w':
		offset_y += 10 * mag;
		break;
	case 'a':
		offset_x -= 10 * mag;
		break;
	case 's':
		offset_y -= 10 * mag;
		break;
	case 'd':
		offset_x += 10 * mag;
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
