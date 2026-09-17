def rank_gf2(matrix):
    if not matrix: return 0
    rows = []
    for r in matrix:
        val = 0
        for bit in r:
            val = (val << 1) | bit
        rows.append(val)
    
    rank = 0
    while rows:
        pivot_row = max(rows)
        if pivot_row == 0:
            break
        rank += 1
        pivot_bit = pivot_row.bit_length() - 1
        
        new_rows = []
        for r in rows:
            if (r >> pivot_bit) & 1:
                r ^= pivot_row
            if r > 0:
                new_rows.append(r)
        rows = new_rows
    return rank

print(rank_gf2([[1,1,0], [0,1,1], [1,0,1]])) # Should be 2, because rows sum to 0 in GF(2)
print(rank_gf2([[1,0,0], [0,1,0], [0,0,1]])) # Should be 3
