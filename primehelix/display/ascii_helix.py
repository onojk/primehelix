import math


def render_ascii_helix(coil, width=64, height=28):
    """
    Double-helix ASCII visualization with:
    - smooth magnitude variation (breathing radius)
    - front/back strand alternation
    - tapered cone shape
    - adaptive connector density
    """

    if coil is None:
        return ["(no helix available)"]

    bit_gap = getattr(coil, "bit_gap", 0)
    balance = getattr(coil, "balance", 0.0)

    # Structural influence
    spread = min(1.0 + balance / 40.0, 3.0)
    compression = max(0.6, 1.0 - bit_gap / 50.0)

    center_x = width // 2
    lines = []

    for y in range(height):
        row = [" "] * width

        # normalized vertical position
        yn = y / max(1, height - 1)

        # phase progression
        t = yn * 2.0 * math.pi / compression

        # smooth "breathing" radius (sin modulation)
        base_r = (width * 0.10) + (width * 0.18) * yn
        breathing = 1.0 + 0.25 * math.sin(2 * t)
        r = base_r * spread * breathing * 0.5

        # two strands (phase offset by pi)
        x1 = int(round(center_x + r * math.cos(t)))
        x2 = int(round(center_x + r * math.cos(t + math.pi)))

        left = min(x1, x2)
        right = max(x1, x2)
        gap = right - left

        # depth alternation (front/back illusion)
        front = (y % 4) < 2

        if 0 <= x1 < width:
            row[x1] = "*" if front else "+"
        if 0 <= x2 < width:
            row[x2] = "+" if front else "*"

        # connector styling
        if 0 <= left < width and 0 <= right < width and gap > 1:
            if y % 2 == 0:
                if gap <= 4:
                    fill = "."
                elif gap <= 12:
                    fill = "~"
                else:
                    fill = "-"
                for x in range(left + 1, right):
                    row[x] = fill

        lines.append("".join(row))

    return lines


def print_ascii_helix(coil):
    lines = render_ascii_helix(coil)

    print()
    print(f"Helix (p={coil.p}, q={coil.q})")
    print(f"balance={coil.balance:.3f}, bit_gap={coil.bit_gap}")
    print()

    for line in lines:
        print(line)
