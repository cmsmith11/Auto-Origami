#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib import collections as mc
import random
import math
import copy

print('starting origami...')

CREASE_SPEED = 2
step = 0.001
EPS = 0.0001
group_colors = ['#73B15B', '#D3974C', '#A467D9', '#E4EE5B', '#EB4343', '#5452E0']
# poly_index = 0

class Circle:
	def __init__(self, i, r, g, x, y):
		self.i = i
		self.r = r # radius
		self.g = g # group
		self.x = x
		self.y = y
		self.links = [] # list of circles an active path away
	def __str__(self):
		return '('+str(self.r)+', '+str(self.g)+')'

class Group:
	def __init__(self, i, x, y):
		self.i = i
		self.x = x
		self.y = y

class Polygon:
	def __init__(self, polypoints):
		self.points = polypoints
		self.unreduced = None # only used for reduced polygons
	def get(self, i):
		for p in self.points:
			if p.i == i:
				return(p)
	def remove(self, i):
		for p in self.points:
			if p.i == i:
				# p.left.right = p.right
				# p.right.left = p.left
				# self.get(p.i).left.right = self.get(p.i).right
				# self.get(p.i).right.left = self.get(p.i).left
				self.points.remove(p)
	def add(self, p):
		self.points.append(p)
		p.left.right = p
		p.right.left = p

class PolyPoint:
	def __init__(self, i, r, g, x, y, left, right):
		self.i = i
		self.r = r
		self.g = g
		self.x = x
		self.y = y
		self.left = left
		self.right = right
		self.theta = None if (left == None or right == None) else find_theta(left, self, right)
	def __str__(self):
		return('Point[' + str(self.i) + ', ' + str(self.r) + ', ' + 
			str(self.g) + ', ' + str(self.x) + ', ' + str(self.y) + ', ' + 
			str(self.left.i) + ', ' + str(self.right.i) + ', ' + str(self.theta) + ']')


# TODOS BY AFTERNOON
# identify all (internal) polygons --> create list <<<<<< DONE
# per polygon
#	create new, h-reduced polygon
# 	step up tiny increments of h, performing certain checks each step
#	each check
#		modify polypoints appropriately (x, y, r) <-- figure out this math for each based on h
#		check between each pair of polypoints in polygon
#		per each polypoint check
#			if true_dist <= group_dist + p1.r + p2.r
#				create two new polygons and recursively split those until polygons degenerate

# does not copy unreduced
def copy_poly(poly):
	polypoints = []
	for p in poly.points:
		polypoints.append(PolyPoint(p.i, p.r, p.g, p.x, p.y, p.left, p.right))
	return(Polygon(polypoints))

# actually checking point lists
def same_polys(poly1, poly2):
	set1 = set([p.i for p in poly1])
	set2 = set([p.i for p in poly2])
	# set1 = set(poly1)
	# set2 = set(poly2)
	return(len(set1.intersection(set2)) == len(set1))

# converts list of circles (all linked properly) to list of polygons
def get_polygons(circles):
	# print('all circles:')
	# for c in circles:
	# 	print([l for l in c.links])

	#  list of list of Circles
	polygon_list = []
	for c in circles:
		for n in c.links:
			poly = get_polygon(c, n, [])
			if poly != None:
				# check if new polygon
				new_poly = True
				for p in polygon_list:
					if same_polys(poly, p):
						new_poly = False
						break
				if new_poly:
					polygon_list.append(poly)
	# list of Polygons
	polygons = []
	# for point_list in polygon_list:
	# 	polypoints = [] # used to make Polygon
	# 	for p in point_list:
	# 		#print('point list:', [point for point in point_list])
	# 		#print('links:', [link for link in p.links])
	# 		ns = []
	# 		for n in point_list:
	# 			if n != p and n in p.links:
	# 				ns.append(n)
	# 		#ns = list(set(point_list).intersection(set(p.links)))
	# 		#print('ns:', ns)
	# 		assert(len(ns) == 2)
	# 		left, right = get_left_right(p, ns[0], ns[1])
	# 		print('left, right:', left, right)
	# 		polypoints.append(PolyPoint(p.i, p.r, p.g, p.x, p.y, left, right))
	# 	polygon = Polygon(polypoints)
	# 	fix_left_rights(polygon)

	# finish left/rights
	for point_list in polygon_list:
		point_list[0].left = point_list[-1]
		for i, p in enumerate(point_list):
			point_list[i].right = point_list[(i + 1) % len(point_list)]
			point_list[i].theta = find_theta(point_list[i].left, point_list[i], point_list[i].right)
		polygon = Polygon(point_list)
		polygons.append(polygon)

	print('all polygons:', len(polygons), polygons)
	for pgon in polygons:
		print('polygon printing...')
		for p in pgon.points:
			print(p)
	return(polygons)


# recursively adds to seen list. Returns None if not internal polygon
def get_polygon(c, n, seen):
	# base case
	for p in seen:
		if c.i == p.i:
			return(seen)
	left = None
	if len(seen) > 0:
		left = seen[-1]
	seen.append(PolyPoint(c.i, c.r, c.g, c.x, c.y, left, None))
	# if c in seen:
	# 	return(seen)
	# seen.append(c)

	ref_p = None
	unit_ns = {} # map from unit neighbor to neighbor
	for nbr in c.links:
		unit_n = norm(c, nbr, 0.2) # only x,y set properly
		unit_ns[unit_n] = nbr
		if nbr == n:
			ref_p = Circle(None, None, None, unit_n.x, unit_n.y)

	while True:
		if (len(unit_ns) < 3):
			return(None)
		rotate_about(c, ref_p, -0.01) # negative for clockwise
		for check_n in unit_ns:
			if unit_ns[check_n] == n:
				continue
			if true_dist(ref_p, check_n) < 0.05:
				return(get_polygon(unit_ns[check_n], c, seen))
		# check if ref_p out of paper??
		if not in_paper(ref_p):
			return(None)

def in_paper(p):
	fudge = 0.1
	return(-fudge < p.x < 1 + fudge and -fudge < p.y < 1 + fudge)

def is_adjacent(p1, p2):
	return(p2.i == p1.left.i or p2.i == p1.right.i)

def unreduced_point(reduced, p):
	for up in reduced.unreduced.points:
		if p.i == up.i:
			return(up)

# return the h-reduced polygon
def hreduce(reduced, h):
	for p in reduced.points:
		a = unreduced_point(reduced, p)
		assert(p.left.i == a.left.i)
		assert(p.right.i == a.right.i)
		#assert(p.theta == find_theta(p.left, p, p.right))
		#print('theta difference?', p.theta, find_theta(a.left, a, a.right))
		p.theta = find_theta(a.left, a, a.right)
		d = h / (math.sin(p.theta / 2))
		new_p = norm(a, avg_p(norm(a, a.left, 10), norm(a, a.right, 10)), d)
		p.x = new_p.x
		p.y = new_p.y
		p.r = a.r - d * math.cos(p.theta / 2)
	return(reduced)

# return two polypoints within distance d, or None
def check_split(polygon, ratio, d):
	#print('checking split!! poly len:', len(polygon.points))
	for p1 in polygon.points:
		for p2 in polygon.points:
			if p1 == p2:
				continue
			if is_adjacent(p1, p2):
				if true_dist(p1, p2) < d:
					print('adjacent close!')
					new_poly = copy_poly(polygon)
					#new_poly.remove(p2.i)
					#left, right = None, None
					if p1.left.i == p2.i:
						new_poly.get(p1.i).left = new_poly.get(p2.i).left
						new_poly.get(p2.left.i).right = new_poly.get(p1.i)
						# p1.left = p2.left
						# p2.left.right = p1
					else:
						new_poly.get(p2.right.i).left = new_poly.get(p1.i)
						new_poly.get(p1.i).right = new_poly.get(p2.i).right
						# p2.right.left = p1
						# p1.right = p2.right
					new_poly.get(p1.i).r = 0
					new_poly.get(p1.i).x = (new_poly.get(p1.i).x + new_poly.get(p2.i).x) / 2
					new_poly.get(p1.i).y = (new_poly.get(p1.i).y + new_poly.get(p2.i).y) / 2
					new_poly.points.remove(new_poly.get(p2.i))
					for np in new_poly.points:
						new_poly.get(np.i).theta = find_theta(np.left, np, np.right)
					# new_p = PolyPoint(p1.i, 0, p1.g, (p1.x + p2.x)/2, (p1.y + p2.y)/2, left, right) # what should radius be??
					# new_poly.add(new_p)
					print('returning new polygon!!', len(polygon.points), len(new_poly.points))
					print('compare orig:')
					for p in polygon.points:
						print(p)
					print('now new:')
					for p in new_poly.points:
						print(p)
					# if len(polygon.points) == 3:
					# 	new_x = 0.5
					# 	new_y = 0.5
					# 	for p in polygon.points:
					# 		polygon.get(p.i).x = new_x
					# 		p.y = new_y
					return([[new_poly]])

			elif true_dist(p1, p2) < unscaled_dist(p1, p2) * ratio + d:
				print('reduced path reached!', true_dist(p1, p2), group_dist(p1.g, p2.g), p1.r, p2.r, d)
				# create two new polygons and recursively split those until polygons degenerate
				new_poly1, new_poly2 = poly_split(polygon, p1, p2)
				return(([new_poly1, new_poly2], [p1, p2]))
				# close.append((p1, p2))
	return(None)

# return two well formatted polygons based on splitting points
def poly_split(polygon, p1, p2):
	#fix_left_rights(polygon)
	print('is well-formatted polygon?')
	for p in polygon.points:
		print(p)
	poly = copy_poly(polygon)
	for i, p in enumerate(poly.points):
		poly.points[i].right = poly.points[(i + 1) % len(poly.points)]
		poly.points[(i + 1) % len(poly.points)].left = poly.points[i]
	print('is well-formatted copied-polygon?')
	for p in poly.points:
		print(p)

	p1_1 = PolyPoint(p1.i, p1.r, p1.g, p1.x, p1.y, poly.get(p1.i).left, None)
	p2_1 = PolyPoint(p2.i, p2.r, p2.g, p2.x, p2.y, p1_1, poly.get(p2.i).right)
	p1_1.right = p2_1
	p1_1.theta = find_theta(p1_1.left, p1, p1_1.right)
	p1_1.left.right = p1_1
	p2_1.right.left = p2_1

	poly.get(p1.i).left = poly.get(p2.i)
	poly.get(p2.i).right = poly.get(p1.i)
	# p1.left = p2
	# p2.right = p1

	poly1_ps = [p1_1]
	curr = p1_1.right
	while curr.i != p1_1.i:
		poly1_ps.append(curr)
		curr = curr.right

	poly2_ps = [poly.get(p2.i)]
	curr = poly.get(p2.i).right
	while curr.i != poly.get(p2.i).i:
		poly2_ps.append(curr)
		curr = curr.right

	ret = [Polygon(poly1_ps), Polygon(poly2_ps)]
	print('check first polygon')
	for p in ret[0].points:
		print(p)
	print('check second polygon')
	for p in ret[1].points:
		print(p)
	return(ret)


def group_dist(g1, g2):
	return group_data[g1][g2] if g1 != g2 else 0

def point_dist(x1, y1, x2, y2):
	return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def true_dist(c1, c2):
	return point_dist(c1.x, c1.y, c2.x, c2.y)
	#return math.sqrt((c2.x - c1.x)**2 + (c2.y - c1.y)**2)

def unscaled_dist(c1, c2):
	return c1.r + c2.r + group_dist(c1.g, c2.g)

# how much the c1 - c2 relation is optimized (number is relative)
def pair_stress(c1, x, y, c2):
	# compute hypothetical distance
	dist = math.sqrt((c2.x - x)**2 + (c2.y - y)**2)
	# compute desired distance
	target = unscaled_dist(c1, c2)
	# compute stress as relative ratio
	return target / dist

# return maximum pairwise stress between all c1 c2 pairs
def total_stress(c, x, y, circles):
	return max([pair_stress(c, x, y, c2) if c != c2 else -1 for c2 in circles]) # is this right??

# clamp val between 0 and 1 inclusive
def clamp(val):
	if val < 0:
		return 0
	elif val > 1:
		return 1
	return val

# return all possible xy locations in neighborhood of check_r
def future_step_options(c):
	check_r = 5 # to check total of (1 + 2*n)**2 spots
	xs = [clamp(c.x + i * step) for i in range(-check_r, check_r + 1)]
	ys = [clamp(c.y + i * step) for i in range(-check_r, check_r + 1)]
	xys = [(x, y) for x in xs for y in ys] # will this work??
	#print(xys)
	return(xys)

# return best future xy location and relative r
def next_xyr(c, circles):
	candidate_xys = future_step_options(c)
	min_stress = None
	champ_xy = (c.x, c.y)
	#print('stresses')
	for x, y in candidate_xys:
		stress = total_stress(c, x, y, circles)
		#print('stress:', stress)
		if min_stress == None or stress < min_stress:
			min_stress = stress
			champ_xy = (x, y)
			#print('new min:', min_stress, champ_xy)
	#print('curr:', (c.x,c.y), 'champ:', champ_xy)
	r = 1 / min_stress
	new_x, new_y = champ_xy
	return((new_x, new_y, r))

# draw circles
def add_circles(ax, circles, ratio):
	for c in circles:
		ax.add_patch(plt.Circle((c.x, c.y), c.r * ratio, color=group_colors[c.g % len(group_colors)]))
		ax.add_patch(plt.Circle((c.x, c.y), 0.01, color='black'))

# animate growing circles. seed_circles' circles are mutated
def animate(ax, seed_circles):
	prev_prev_xyrs = []
	prev_xyrs = []
	min_r = None
	while True:
		# for each circle c1, find step of least stress
		next_xyrs = [next_xyr(c, seed_circles) for c in seed_circles]
		# adjust circles' locations
		min_r = min([xyr[2] for xyr in next_xyrs])
		for i, c in enumerate(seed_circles):
			c.x, c.y, _ = next_xyrs[i]

		# display current time step with proper scaling
		add_circles(ax, seed_circles, min_r)

		# break if stopped
		if next_xyrs == prev_xyrs or next_xyrs == prev_prev_xyrs:
			print(min_r)
			break
		prev_prev_xyrs = prev_xyrs
		prev_xyrs = next_xyrs

		# pause and move to next time step
		plt.pause(0.001)
		ax.clear()
	print('out of while loop')

	plt.pause(0.5)
	draw_creases(ax, seed_circles, min_r)
	plt.pause(0.001) # neccessary for some reason

	# plt.show()

def draw_creases(ax, circles, ratio):
	print('drawing creases!')
	axials = []
	ridges = []
	for c1 in circles:
		for c2 in circles:
			if c1 != c2 and (true_dist(c1, c2) - unscaled_dist(c1, c2) * ratio < 0.05): # can tweak this constant
				axials.append([(c1.x, c1.y), (c2.x, c2.y)])
				c1.links.append(c2)
				c2.links.append(c1)
				print('linking!')
	# link side circles ?
	buf = 0.02
	top = sorted([c for c in circles if c.y > 1 - buf], key=lambda c : c.x)
	right = sorted([c for c in circles if c.x > 1 - buf], key=lambda c : c.y)
	bottom = sorted([c for c in circles if c.y < buf], key=lambda c : c.x)
	left = sorted([c for c in circles if c.x < buf], key=lambda c : c.y)

	all_sides = [top, right, bottom, left]
	print('allsides:', all_sides)
	for side in all_sides:
		for i in range(len(side) - 1):
			if side[i + 1] not in side[i].links:
				print('adding side1!!')
				side[i].links.append(side[i + 1])
			if side[i] not in side[i + 1].links:
				side[i + 1].links.append(side[i])
				print('adding side2!!')





	print(len(axials))
	lc = mc.LineCollection(axials, color='red')
	ax.add_collection(lc)

	plt.pause(0.5)

	# get polygons
	polygons = get_polygons(circles)
	#print('polygons??', len(polygons), polygons)
	mountains = []

	# per polygon
	for poly in polygons:
		print('new polygon!!')
		for p in poly.points:
			print(p)
		# create new, h-reduced polygon
		h = 0
		reduced = copy_poly(poly)
		reduced.unreduced = poly
		# step up tiny increments of h, performing certain checks each step
		poly_looping = True
		while poly_looping:
			h += 0.001 * CREASE_SPEED
			# modify polypoints appropriately (x, y, r)
			reduced = hreduce(reduced, h)
			# check if reduced points close to any other
			split_check = check_split(reduced, ratio, 0.006)
			if split_check != None:
				if len(split_check) > 1:
					# contains polys and split points
					p1, p2 = split_check[1]
					mountains.append([(p1.x, p1.y), (p2.x, p2.y)])
				for new_poly in split_check[0]:
					if len(new_poly.points) >= 3: # ignore degenerate polygons
						polygons.append(new_poly)
					# elif len(new_poly.points) == 2:
					# 	p1, p2 = new_poly.points
					# 	ridges.append([(p1.x, p1.y), (p2.x, p2.y)])
				if len(reduced.points) == 3:
					p1, p2, p3 = reduced.points
					new_x = (p1.x + p2.x + p3.x) / 3
					new_y = (p1.y + p2.y + p3.y) / 3
					for p in reduced.points:
						p.x, p.y = new_x, new_y
				poly_looping = False
				#break

			for rp in reduced.points:
				urp = unreduced_point(reduced, rp)
				ridges.append([(rp.x, rp.y), (urp.x, urp.y)])

			lc2 = mc.LineCollection(ridges, color='blue')
			lc3 = mc.LineCollection(mountains, color='red')
			add_circles(ax, circles, ratio)
			ax.add_collection(lc3)
			ax.add_collection(lc2)
			ax.add_collection(lc)
			plt.pause(0.0001)
			ax.clear()

		# add lines to draw
		# for rp in reduced.points:
		# 	# add line from unreduced to reduced
		# 	urp = unreduced_point(reduced, rp)
		# 	ridges.append([(rp.x, rp.y), (urp.x, urp.y)])


	#		per each polypoint check
	#			if true_dist <= group_dist + p1.r + p2.r
	#				create two new polygons and recursively split those until polygons degenerate

	add_circles(ax, circles, ratio)
	ax.add_collection(lc3)
	ax.add_collection(lc2)
	ax.add_collection(lc)
	# lc2 = mc.LineCollection(ridges)
	# ax.add_collection(lc2)

	#add_circles(ax, circles, ratio)
	print('drew creases!')

	return(circles)
	# plt.show()

# p1, vert, p3 are three neighboring Circles
def find_theta(p1, vert, p3):
	# law of cosines
	a = true_dist(p1, vert)
	b = true_dist(p3, vert)
	c = true_dist(p1, p3)
	theta = math.acos((a**2 + b**2 - c**2) / (2 * a * b))
	return(theta)

# takes two Circles
def avg_p(p1, p2):
	return(Circle(None, None, None, (p1.x + p2.x)/2, (p1.y+p2.y)/2))

# takes two Circles
def norm(tail, tip, scale=1):
	dx = tip.x - tail.x
	dy = tip.y - tail.y
	d = true_dist(tail, tip)
	return(Circle(None, None, None, tail.x + dx/d*scale, tail.y + dy/d*scale))

# both Circles, theta is in radians, rotates counter-clockwise
def rotate_about(center, p, theta):
	new_x = (p.x - center.x) * math.cos(theta) - (p.y - center.y) * math.sin(theta) + center.x
	new_y = (p.x - center.x) * math.sin(theta) + (p.y - center.y) * math.cos(theta) + center.y
	p.x, p.y = new_x, new_y

def mag(vec):
	sqrSum = 0
	for i in vec:
		sqrSum += i**2
	return(math.sqrt(sqrSum))

# takes 2 3d vectors
def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
    return c

def get_left_right(c, n1, n2):
	a = (n1.x - c.x, n1.y - c.y, 0)
	b = (n2.x - c.x, n2.y - c.y, 0)
	cp = cross(a, b)
	# p1 = PolyPoint(n1.i, n1.r, n1.g, n1.x, n1.y, None, None)
	# p2 = PolyPoint()
	if mag(cp) > 0:
		return((n1, n2)) # n1 is left
	else:
		return((n2, n1)) # n2 is left

def fix_left_rights(polygon):
	for p in polygon.points:
		p.left = polygon.get(p.left.i)
		p.right = polygon.get(p.right.i)



# lines = [[(0, 1), (1, 1)], [(2, 3), (3, 3)], [(1, 2), (1, 3)]]
# c = np.array([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)])

# lc = mc.LineCollection(lines, colors=c, linewidths=2)
# fig, ax = pl.subplots()
# ax.add_collection(lc)



# read .stk format and put into point_info: {i: [x, y, [n1, n2, ...]], ...}
def stk_to_dict(filename):
	point_info = {}
	f = open(filename, 'r')
	for i, line in enumerate(f):
		line = line[:-1]
		print(line)
		parts = line.split(' ')
		point_info[i] = [float(parts[0]), float(parts[1])]
		ns = parts[2].split(',')
		point_info[i].append([int(n) for n in ns])
	f.close()
	print(point_info)
	return(point_info)

# convert stk dict to 3 data structures: circle_data, group_data and g_locs
# def dict_to_structures(point_info):
# 	circle_data = []
# 	group_data = {}
# 	g_locs = {}

# 	# if at least one of a point's neighbors is a leaf, it is a group

# 	return((circle_data, group_data, g_locs))

# recursively find distance from g1 to g2 through links
def find_group_dist(g1, g2, from_p, point_info):
	# for each g1 neighbor
	x1, y1 = point_info[g1][0], point_info[g1][1]
	for n in point_info[g1][2]:
		if n == from_p:
			continue
		if n == g2:
			x2, y2 = point_info[n][0], point_info[n][1]
			return(point_dist(x1, y1, x2, y2))
		x2, y2 = point_info[n][0], point_info[n][1]
		ret = find_group_dist(n, g2, g1, point_info)
		if ret != None:
			return(point_dist(x1, y1, x2, y2) + ret)
	return(None)

def dict_to_structures(point_info):
	groups = []
	for p in point_info:
		is_group = False
		for n in point_info[p][2]:
			if len(point_info[n][2]) == 1:
				is_group = True
				break
		if is_group:
			groups.append(p)
	print('groups', groups)

	g_locs = {}
	for i, g in enumerate(groups):
		g_locs[i] = (point_info[g][0], point_info[g][1])

	circle_data = []
	for i, g in enumerate(groups):
		x1, y1 = point_info[g][0], point_info[g][1]
		circle_group = [point_dist(x1, y1, point_info[n][0], point_info[n][1]) for n in point_info[g][2] if n not in groups]
		circle_data.append(circle_group)

	group_data = {}
	for i, g1 in enumerate(groups):
		group_data[i] = {}
		for j, g2 in enumerate(groups):
			if i != j:
				group_data[i][j] = find_group_dist(g1, g2, None, point_info)

	# return((circle_data, group_data, g_locs))
	return((circle_data, group_data, groups))


# point_info = stk_to_dict('temp.stk')
# point_info = stk_to_dict('insect.stk')
# point_info = stk_to_dict('10-flap-leaf.stk')
# point_info = stk_to_dict('1-7.stk')
# point_info = stk_to_dict('8-equal-flaps.stk')
# point_info = stk_to_dict('man.stk')
# point_info = stk_to_dict('hand.stk')
# point_info = stk_to_dict('dragon.stk')
# point_info = stk_to_dict('long-tail-rat.stk')
point_info = stk_to_dict('3-short-front-3-long-back.stk')
#circle_data, group_data, g_locs = dict_to_structures(point_info)
circle_data, group_data, groups = dict_to_structures(point_info)

print('circle_data:', circle_data)
print('group_data:', group_data)
#print('g_locs:', g_locs)


# circle_data = [
# 	[3], # head
# 	[1, 1, 1, 1, 1], # hand
# 	[1, 1, 1, 1, 1], # hand
# 	[8, 8] # legs
# ]

# group_data = {
# 	0: {3: 5, 2: 7, 1: 7},
# 	1: {0: 7, 3: 12, 2: 14},
# 	2: {0: 7, 3: 12, 1: 14},
# 	3: {0: 5, 2: 12, 1: 12}
# }

# g_locs = {
# 	0: (0.99, 0.99),
# 	1: (0.01, 0.99),
# 	2: (0.99, 0.01),
# 	3: (0.01, 0.01)
# }

# circle_data = [
# 	[4, 4, 4],
# 	[2, 3, 2]
# ]

# group_data = {
# 	0: {1: 1},
# 	1: {0: 1}
# }

# circle_data = [
# 	[1 for i in range(7)]
# ]

# group_data = {}

# 0.) recieve data in p(x, y, [neighbors]) form, convert
# 1.) normalize lengths
# 2.) 

# that's syntatically correct?
# so i goes 0, 1, 2, 3 yes?
# g is circle data so it goes
# [circle (4, 0)(4,0)(4,0)
# cirlce (2,1) (3, 1) (2,1)]




# circles = []
# for p in point_info:
# 	if p not in groups:
# 		x1, y1 = point_info[p][0], point_info[p][1]
# 		n = point_info[p][2][0]
# 		x2, y2 = point_info[n][0], point_info[n][1]
# 		circles.append(Circle(point_dist(x1, y1, x2, y2), Group(groups.index(n), -1, -1), x1, y1))

#circles = [Circle(r, Group(i, g_locs[i][0], g_locs[i][1])) for i, g in enumerate(circle_data) for r in g]

fig, ax = plt.subplots()
fig.set_size_inches(6, 6)

num_trials = 1
max_scale = None
champ_circles = None
champ_seed = None
for i in range(num_trials):
	# circles = [Circle(r, i) for i, g in enumerate(circle_data) for r in g]
	#circles = [Circle(r, Group(i, g_locs[i][0], g_locs[i][1])) for i, g in enumerate(circle_data) for r in g]
	circles = []
	for p in point_info:
		if p not in groups:
			x1, y1 = point_info[p][0], point_info[p][1]
			n = point_info[p][2][0]
			x2, y2 = point_info[n][0], point_info[n][1]
			#circles.append(Circle(point_dist(x1, y1, x2, y2), Group(groups.index(n), -1, -1), x1, y1))
			circles.append(Circle(p, point_dist(x1, y1, x2, y2), groups.index(n), x1, y1))
	seed = copy.deepcopy(circles)

	prev_prev_xyrs = []
	prev_xyrs = []
	# trial starting
	while True:
		next_xyrs = [next_xyr(c, circles) for c in circles]

		# adjust circle locations
		for i, c in enumerate(circles):
			c.x, c.y, _ = next_xyrs[i]

		# record min scale
		min_r = min([xyr[2] for xyr in next_xyrs])

		# break if stopped
		if next_xyrs == prev_xyrs or next_xyrs == prev_prev_xyrs:
			break
		prev_prev_xyrs = prev_xyrs
		prev_xyrs = next_xyrs

	if True:#if max_scale == None or min_r > max_scale:
		max_scale = min_r
		champ_circles = circles
		champ_seed = seed

#add_circles(ax, champ_circles, max_scale)
print(max_scale)
animate(ax, champ_seed)
print('finished animation')
plt.show()
