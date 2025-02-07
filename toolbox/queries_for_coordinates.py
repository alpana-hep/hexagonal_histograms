#!/usr/bin/env python2

query_line1 = [
	[[8, 18], [1, 0]],
	[[8, 18, 17], [2, 5, 1]],
	[[18, 27, 17], [4, 0, 2]],
	[[17, 26, 27], [3, 1, 5]],
	[[26, 27, 37], [2, 4, 0]],
	[[26, 36, 37], [3, 1, 5]],
	[[36, 37, 48], [2, 4, 0]],
	[[36, 47, 48], [3, 1, 5]],
	[[47, 48, 60], [2, 4, 0]],
	[[47, 59, 60], [3, 1, 5]],
	[[59, 60, 75], [2, 4, 0]],
	[[59, 74, 75], [3, 1, 5]],
	[[74, 75, 89], [2, 4, 0]],
	[[74, 88, 89], [3, 1, 5]],
	[[88, 89, 104], [2, 4, 0]],
	[[88, 103, 104], [3, 1, 5]],
]

query_line2 = [
	[[81, 96], [4, 0]],
	[[81, 96, 97], [3, 1, 5]],
	[[81, 82, 97], [2, 4, 0]],
	[[82, 97, 98], [3, 1, 5]],
	[[82, 83, 98], [2, 4, 0]],
	[[83, 98, 99], [3, 1, 5]],
	[[83, 84, 99], [2, 4, 0]],
	[[84, 99, 100], [3, 1, 5]],
	[[84, 85, 100], [2, 4, 0]],
	[[85, 100, 101], [3, 1, 5]],
	[[85, 86, 101], [2, 4, 0]],
	[[86, 101, 102], [3, 1, 5]],
	[[86, 87, 102], [2, 4, 0]],
	[[87, 102, 103], [3, 1, 5]],
	[[87, 88, 103], [2, 4, 0]],
	[[88, 103, 104], [3, 1, 5]],
]

query_line3 = [
	[[88, 103, 104], [3, 1, 5]],
	[[103, 104, 119], [2, 4, 0]],
	[[104, 119, 120], [3, 1, 5]],
	[[119, 120, 134], [2, 4, 0]],
	[[120, 134, 135], [3, 1, 5]],
	[[134, 135, 149], [2, 4, 0]],
	[[135, 149, 150], [3, 1, 5]],
	[[149, 150, 164], [2, 4, 0]],
	[[150, 164, 165], [3, 1, 5]],
	[[164, 165, 176], [2, 4, 0]],
	[[165, 176, 177], [3, 1, 5]],
	[[176, 177, 187], [2, 4, 0]],
	[[177, 187, 188], [3, 1, 5]],
	[[187, 188, 197], [2, 4, 0]],
	[[188, 197, 198], [3, 1, 4]],
	[[197, 198], [2, 3]],
]

query_line4 = [
	[[29, 40], [4, 0]],
	[[29, 40, 41], [3, 1, 5]],
	[[29, 30, 41], [2, 4, 0]],
	[[30, 41, 42], [3, 1, 5]],
	[[30, 31, 42], [2, 4, 0]],
	[[31, 42, 43], [3, 1, 5]],
	[[31, 32, 43], [2, 4, 0]],
	[[32, 43, 44], [3, 1, 5]],
	[[32, 33, 44], [2, 4, 0]],
	[[33, 44, 45], [3, 1, 5]],
	[[33, 34, 45], [2, 4, 0]],
	[[34, 45, 46], [3, 1, 5]],
	[[34, 35, 46], [2, 4, 0]],
	[[35, 46, 47], [3, 1, 5]],
	[[35, 36, 47], [2, 4, 0]],
	[[36, 47, 48], [3, 1, 5]],
]

query_line5 = [
	[[51, 65], [1, 0]],
	[[51, 64, 65], [2, 1, 4]],
	[[64, 65, 79], [2, 3, 0]],
	[[64, 78, 79], [3, 1, 5]],
	[[78, 79, 93], [2, 4, 0]],
	[[78, 92, 93], [3, 1, 5]],
	[[92, 93, 108], [2, 4, 0]],
	[[92, 107, 108], [3, 1, 5]],
	[[107, 108, 123], [2, 4, 0]],
	[[107, 122, 123], [3, 1, 5]],
	[[122, 123, 137], [2, 4, 0]],
	[[122, 136, 137], [3, 1, 5]],
	[[136, 137, 151], [2, 4, 0]],
	[[136, 150, 151], [3, 1, 5]],
	[[150, 151, 165], [2, 4, 0]],
	[[150, 164, 165], [3, 1, 5]],
]

query_line6 = [
	[[84, 99, 100], [3, 1, 5]],
	[[99, 100, 115], [2, 4, 0]],
	[[100, 115, 116], [3, 1, 5]],
	[[115, 116, 130], [2, 4, 0]],
	[[116, 130, 131], [3, 1, 5]],
	[[130, 131, 145], [2, 4, 0]],
	[[131, 145, 146], [3, 1, 5]],
	[[145, 146, 159], [2, 4, 0]],
	[[146, 159, 160], [3, 1, 5]],
	[[159, 160, 172], [2, 4, 0]],
	[[160, 172, 173], [3, 1, 5]],
	[[172, 173, 183], [2, 4, 0]],
	[[173, 183, 184], [3, 1, 5]],
	[[183, 184, 193], [2, 4, 0]],
	[[184, 193, 194], [3, 1, 4]],
	[[193, 194], [2, 3]],
]
