#include <GL/glut.h>
#include <stdio.h>
#include <math.h>
#include <queue>
#include <limits>
#include <random>
#include "biarc_approx.h"
#include "aabb.h"

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
bool isDrawAABB = true;
bool isDrawHierarchy = false;
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

void draw_line(const Point p, const Point tan){
	Point pt;
	glBegin(GL_LINE_STRIP);
	pt[0] = p[0] - tan[0] * width * 10.0;
	pt[1] = p[1] - tan[1] * height * 10.0;
	glVertex2f(pt[0], pt[1]);
	pt[0] = p[0] + tan[0] * width * 10.0;
	pt[1] = p[1] + tan[1] * height * 10.0;
	glVertex2f(pt[0], pt[1]);
	glEnd();
}

void build_hierarchy(std::shared_ptr<Hierarchy> h, int power){
	std::vector<CubicBezierCurve> curves;
	CubicBezierCurve &seg = h->curve;

	// divide seg at inflection point
	std::vector<CubicBezierCurve> segs;
	subdivide(&seg, segs, 1);

	auto leftH = std::make_shared<Hierarchy>();
	auto rightH = std::make_shared<Hierarchy>();
	leftH->curve = segs[0];
	rightH->curve = segs[1];

	if (power == 0){
		Point inflect;
		Arc arc1, arc2;

		get_biarc_inflect(&seg, inflect);
		to_biarc(&seg, inflect, &arc1, &arc2);

		Arc line1, line2;
		copy_point(seg.control_pts[0], line1.center);
		copy_point(seg.control_pts[3], line2.center);
		line1.radius = NAN;
		line2.radius = NAN;
		line1.begin = inflect[0];
		line1.end = inflect[1];
		line2.begin = inflect[0];
		line2.end = inflect[1];
		
		REAL error1 = get_AABB(segs[0], arc1, leftH->box, isDrawBiarcs && !power);
		REAL error2 = get_AABB(segs[1], arc2, rightH->box, isDrawBiarcs && !power);

		AABB box1, box2;
		REAL box_e1 = get_AABB(segs[0], line1, box1, isDrawBiarcs && !power);
		REAL box_e2 = get_AABB(segs[1], line2, box2, isDrawBiarcs && !power);

		if (box_e1 < error1){
			arc1 = line1;
			leftH->box = box1;
			error1 = box_e1;
		}
		if (box_e2 < error2){
			arc2 = line2;
			rightH->box = box2;
			error2 = box_e2;
		}

		leftH->box.x[0] -= error1;
		leftH->box.x[1] += error1;
		leftH->box.y[0] -= error1;
		leftH->box.y[1] += error1;
		rightH->box.x[0] -= error2;
		rightH->box.x[1] += error2;
		rightH->box.y[0] -= error2;
		rightH->box.y[1] += error2;
		leftH->arc = std::make_shared<Arc>(arc1);
		rightH->arc = std::make_shared<Arc>(arc2);
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

REAL sample_points_distance(CubicBezierCurve c1, CubicBezierCurve c2, int num_samples){
	std::vector<REAL> pts1x, pts1y;
	std::vector<REAL> pts2x, pts2y;
	for (int i = 0; i <= num_samples; i++){
		Point pt;
		const REAL t = (REAL)i / (REAL)num_samples;
		evaluate(&c1, t, pt);
		pts1x.push_back(pt[0]);
		pts1y.push_back(pt[1]);
	}
	for (int i = 0; i <= num_samples; i++){
		Point pt;
		const REAL t = (REAL)i / (REAL)num_samples;
		evaluate(&c2, t, pt);
		pts2x.push_back(pt[0]);
		pts2y.push_back(pt[1]);
	}

	REAL min_dist = std::numeric_limits<REAL>::max();
	for (int i = 0; i < pts1x.size(); i++){
		for (int j = 0; j < pts2x.size(); j++){
			Point pt1, pt2;
			pt1[0] = pts1x[i];
			pt1[1] = pts1y[i];
			pt2[0] = pts2x[j];
			pt2[1] = pts2y[j];
			const REAL dist = distance(pt1, pt2);
			if (dist < min_dist){
				min_dist = dist;
			}
		}
	}

	return min_dist;
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

typedef std::pair<REAL, std::pair<std::shared_ptr<Hierarchy>, std::shared_ptr<Hierarchy>>> min_pair;
void draw_minimum_distance(std::shared_ptr<Hierarchy> tree1, std::shared_ptr<Hierarchy> tree2){
	std::priority_queue<min_pair, std::vector<min_pair>, std::greater<min_pair>> q;
	
	// Use bounding box for bound computation, use biarc for final computation
	REAL lower_bound = distance(tree1->box, tree2->box);
	REAL upper_bound = sample_points_distance(tree1->curve, tree2->curve, NUM_SAMPLES);
	CubicBezierCurve& bound_curve1 = tree1->curve, &bound_curve2 = tree2->curve;

	q.push(std::make_pair(lower_bound, std::make_pair(tree1, tree2)));
	while (!q.empty()){
		REAL curr_bound = q.top().first;
		if (curr_bound > upper_bound)
			break;
		lower_bound = curr_bound;
		auto tree1 = q.top().second.first;
		auto tree2 = q.top().second.second;
		q.pop();

		REAL local_bound = sample_points_distance(tree1->curve, tree2->curve, NUM_SAMPLES);
		if (upper_bound > local_bound) {
			upper_bound = local_bound;
			bound_curve1 = tree1->curve;
			bound_curve2 = tree2->curve;
		}

		if (tree1->left == nullptr && tree2->left == nullptr){
			// Both BVH reached leaf node
			// Set arc distance as upper bound, ignore biarc approximation error
			REAL local_distance = distance(tree1->arc, tree2->arc);
			if (local_distance < upper_bound){
				if (local_distance < lower_bound){
					local_distance = lower_bound;
				}
				upper_bound = local_distance;
				bound_curve1 = tree1->curve;
				bound_curve2 = tree2->curve;
			}
		}

		if (upper_bound < lower_bound){
			lower_bound = upper_bound;
			break;
		}

		// Add child nodes to priority queue
		if (tree1->left != nullptr && (volume(tree1->box) < volume(tree2->box) || tree2->left == nullptr)){
			auto l_lower_bound = distance(tree1->left->box, tree2->box);
			if (l_lower_bound < upper_bound){
				q.push(std::make_pair(l_lower_bound, std::make_pair(tree1->left, tree2)));
			}
			auto r_lower_bound = distance(tree1->right->box, tree2->box);
			if (r_lower_bound < upper_bound){
				q.push(std::make_pair(r_lower_bound, std::make_pair(tree1->right, tree2)));
			}
		}
		else if (tree2->left != nullptr && (volume(tree1->box) > volume(tree2->box) || tree1->left == nullptr)){
			auto l_lower_bound = distance(tree1->box, tree2->left->box);
			if (l_lower_bound < upper_bound){
				q.push(std::make_pair(l_lower_bound, std::make_pair(tree1, tree2->left)));
			}
			auto r_lower_bound = distance(tree1->box, tree2->right->box);
			if (r_lower_bound < upper_bound){
				q.push(std::make_pair(r_lower_bound, std::make_pair(tree1, tree2->right)));
			}
		}
	}
	
	std::string distance = "Distance: " + std::to_string((upper_bound + lower_bound) / 2);
	std::string error = "Error: " + std::to_string((upper_bound - lower_bound) / 2.0);
	std::string arc_counts = "Num of Arcs: " + std::to_string(int(std::pow(2.0, subdivision_power)));
	glLineWidth(10.0);
	glColor3ub(255, 0, 0);
	draw_curve(bound_curve1);
	draw_curve(bound_curve2);
	glLineWidth(1.0);
	draw_text(arc_counts);
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

	auto root1 = std::make_shared<Hierarchy>();
	auto root2 = std::make_shared<Hierarchy>();

	root1->curve = curve1;
	build_hierarchy(root1, subdivision_power);
	root2->curve = curve2;
	build_hierarchy(root2, subdivision_power);

	draw_minimum_distance(root1, root2);
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
		// isDrawBiarcs ^= true;
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
