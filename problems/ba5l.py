import sys


file = open("rosalind_ba5l.txt")
s, t = [line.strip() for line in file.readlines()]

score_mat = [  # substitution penalty
	[4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2],
	[0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
	[-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3],
	[-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2],
	[-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3],
	[0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3],
	[-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2],
	[-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1],
	[-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2],
	[-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1],
	[-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1],
	[-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2],
	[-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3],
	[-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1],
	[-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2],
	[1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2],
	[0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2],
	[0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1],
	[-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2],
	[-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7]
]

score_mat_map = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
	'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}


def hirschberg(a, b, la, ra, lb, rb, score_mat, score_mat_map, gamma):
	"""
		e a1 a2 a3 ... am
	e
	b1
	b2
	b3
	...
	bn
	"""
	assert ra - la > 0
	assert rb - lb > 0
	if ra - la == 1 and rb - lb == 1:
		return a[la], b[lb]
	elif ra - la == 1:
		max_subst_cost = -sys.maxsize
		max_subst_b_idx = None
		for j in range(lb, rb):
			tmp = score_mat[score_mat_map[a[la]]][score_mat_map[b[j]]]
			if tmp > max_subst_cost:
				max_subst_cost = tmp
				max_subst_b_idx = j
		if max_subst_cost < gamma:
			return a[la] + '-' * (rb - lb), '-' + b[lb:rb]
		else:
			return '-' * (max_subst_b_idx - lb) + a[la] + '-' * (rb - max_subst_b_idx - 1), b[lb:rb]
	elif rb - lb == 1:
		max_subst_cost = -sys.maxsize
		max_subst_a_idx = None
		for i in range(la, ra):
			tmp = score_mat[score_mat_map[a[i]]][score_mat_map[b[lb]]]
			if tmp > max_subst_cost:
				max_subst_cost = tmp
				max_subst_a_idx = i
		if max_subst_cost < gamma:
			return '-' + a[la:ra], b[lb] + '-' * (ra - la)
		else:
			return a[la:ra], '-' * (max_subst_a_idx - la) + b[lb] + '-' * (ra - max_subst_a_idx - 1)
	else:
		mid_b = (lb + rb) // 2  # Horizontal mid line.
		prefix_dists = [j * gamma for j in range(ra - la + 1)]
		for i in range(mid_b - lb):
			old_dist = prefix_dists[0]
			prefix_dists[0] = (i + 1) * gamma
			for j in range(ra - la):
				subst_cost = score_mat[score_mat_map[a[la + j]]][score_mat_map[b[lb + i]]]
				tmp = max(old_dist + subst_cost, prefix_dists[j] + gamma, prefix_dists[j + 1] + gamma)
				old_dist = prefix_dists[j + 1]
				prefix_dists[j + 1] = tmp
		suffix_dists = [j * gamma for j in range(ra - la + 1)][::-1]
		for i in reversed(range(rb - mid_b)):
			old_dist = suffix_dists[-1]
			suffix_dists[-1] = (rb - mid_b - i) * gamma
			for j in reversed(range(ra - la)):
				subst_cost = score_mat[score_mat_map[a[la + j]]][score_mat_map[b[mid_b + i]]]
				tmp = max(old_dist + subst_cost, suffix_dists[j] + gamma, suffix_dists[j + 1] + gamma)
				old_dist = suffix_dists[j]
				suffix_dists[j] = tmp
		mid_a = None
		max_a_intersection_val = -sys.maxsize
		assert len(prefix_dists) == len(suffix_dists)
		for idx in range(len(prefix_dists)):
			tmp = prefix_dists[idx] + suffix_dists[idx]
			if prefix_dists[idx] + suffix_dists[idx] > max_a_intersection_val:
				max_a_intersection_val = tmp
				mid_a = la + idx
		if mid_a > la:
			a_left, b_left = hirschberg(s, t, la, mid_a, lb, mid_b, score_mat, score_mat_map, gamma)
		else:
			a_left = '-' * (mid_b - lb)
			b_left = b[lb:mid_b]
		if ra > mid_a:
			a_right, b_right = hirschberg(s, t, mid_a, ra, mid_b, rb, score_mat, score_mat_map, gamma)
		else:
			a_right = '-' * (rb - mid_b)
			b_right = b[mid_b:rb]
		return a_left + a_right, b_left + b_right


s_align, t_align = hirschberg(s, t, 0, len(s), 0, len(t), score_mat, score_mat_map, -5)
assert len(s_align) == len(t_align)
cost = 0
for idx in range(len(s_align)):
	if s_align[idx] == '-' or t_align[idx] == '-':
		cost -= 5
	else:
		cost += score_mat[score_mat_map[s_align[idx]]][score_mat_map[t_align[idx]]]
print(cost)
print(s_align)
print(t_align)