with open("rosalind_hamm.txt") as file:
	lines = file.readlines()

s = lines[0].strip()
t = lines[1].strip()

assert len(s) == len(t)

distance = 0
for i in range(len(s)):
	if s[i] != t[i]:
		distance += 1

print(distance)