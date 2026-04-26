from manim import *

SLOW = 4  # 4x slower animation timing


class SemiprimeLopsidedExplainer(Scene):
    def construct(self):
        self.camera.background_color = BLACK

        def slow_play(*anims, run_time=1):
            self.play(*anims, run_time=run_time * SLOW)

        def slow_wait(t=1):
            self.wait(t * SLOW)

        # ------------------------------------------------------------------
        # Title
        # ------------------------------------------------------------------
        title = Text("Why are most semiprimes lopsided?", font_size=42)
        subtitle = Text("Small primes divide many numbers", font_size=30)
        subtitle.next_to(title, DOWN)

        slow_play(Write(title), FadeIn(subtitle, shift=DOWN))
        slow_wait(1)
        slow_play(FadeOut(title), FadeOut(subtitle))

        # ------------------------------------------------------------------
        # Number line / divisibility
        # ------------------------------------------------------------------
        line = NumberLine(
            x_range=[0, 40, 5],
            length=11,
            include_numbers=True,
            font_size=22,
        ).shift(UP * 0.25)

        label = Text("Numbers from 1 to 40", font_size=30).next_to(line, UP, buff=0.7)

        slow_play(Create(line), Write(label))
        slow_wait(0.5)

        dots_2 = VGroup(*[Dot(line.n2p(i), radius=0.07) for i in range(2, 41, 2)])
        dots_3 = VGroup(*[Dot(line.n2p(i), radius=0.06) for i in range(3, 41, 3)])
        dots_5 = VGroup(*[Dot(line.n2p(i), radius=0.055) for i in range(5, 41, 5)])

        text_small = Text(
            "A prime “hits” a number when it divides it evenly",
            font_size=30,
        ).to_edge(UP)

        slow_play(ReplacementTransform(label, text_small))
        slow_wait(0.5)

        caption_2 = Text("2 divides about 1/2 of all numbers", font_size=30).to_edge(DOWN)
        caption_3 = Text("3 divides about 1/3 of all numbers", font_size=30).to_edge(DOWN)
        caption_5 = Text("5 divides about 1/5 of all numbers", font_size=30).to_edge(DOWN)

        slow_play(FadeIn(dots_2), Write(caption_2))
        slow_wait(0.75)
        slow_play(ReplacementTransform(caption_2, caption_3), FadeIn(dots_3))
        slow_wait(0.75)
        slow_play(ReplacementTransform(caption_3, caption_5), FadeIn(dots_5))
        slow_wait(0.75)

        explanation = Text(
            "Small primes divide many numbers.\nLarge primes divide very few.",
            font_size=34,
            line_spacing=0.9,
        ).to_edge(DOWN)

        slow_play(ReplacementTransform(caption_5, explanation))
        slow_wait(1.25)

        slow_play(
            FadeOut(dots_2),
            FadeOut(dots_3),
            FadeOut(dots_5),
            FadeOut(explanation),
            FadeOut(line),
            FadeOut(text_small),
        )

        # ------------------------------------------------------------------
        # Semiprime explanation
        # ------------------------------------------------------------------
        semiprime_title = Text("A semiprime has exactly two prime factors", font_size=36)
        semiprime_expr = Text("semiprime = prime × prime", font_size=38)
        semiprime_expr.next_to(semiprime_title, DOWN, buff=0.6)

        slow_play(Write(semiprime_title), FadeIn(semiprime_expr, shift=DOWN))
        slow_wait(1.25)
        slow_play(FadeOut(semiprime_title), FadeOut(semiprime_expr))

        # ------------------------------------------------------------------
        # Common vs rare comparison
        # ------------------------------------------------------------------
        left_title = Text("Common", font_size=34)
        left_expr = Text("3 × 1,000,003", font_size=44)
        left_note = Text("small × large", font_size=28)

        left_group = VGroup(left_title, left_expr, left_note)
        left_group.arrange(DOWN, buff=0.35)
        left_group.move_to(LEFT * 3.35 + UP * 1.2)

        right_title = Text("Rare", font_size=34)
        right_expr = Text("1009 × 1013", font_size=44)
        right_note = Text("large × large", font_size=28)

        right_group = VGroup(right_title, right_expr, right_note)
        right_group.arrange(DOWN, buff=0.35)
        right_group.move_to(RIGHT * 3.35 + UP * 1.2)

        slow_play(FadeIn(left_group, shift=UP), FadeIn(right_group, shift=UP))
        slow_wait(0.5)

        arrow = Arrow(
            start=left_expr.get_right() + RIGHT * 0.35,
            end=right_expr.get_left() + LEFT * 0.35,
            buff=0.15,
            stroke_width=6,
            max_tip_length_to_length_ratio=0.18,
        )

        comparison = Text(
            "Small × large has many more chances.\nLarge × large needs two rare events at once.",
            font_size=30,
            line_spacing=0.8,
        ).to_edge(DOWN)

        slow_play(GrowArrow(arrow))
        slow_play(Write(comparison))
        slow_wait(1.5)

        slow_play(
            FadeOut(left_group),
            FadeOut(right_group),
            FadeOut(arrow),
            FadeOut(comparison),
        )

        # ------------------------------------------------------------------
        # Divisor density summary
        # ------------------------------------------------------------------
        density = Text(
            "The key idea:\nsemiprimes are shaped by divisibility,\nnot by evenly pairing primes.",
            font_size=36,
            line_spacing=0.9,
        )

        slow_play(Write(density))
        slow_wait(1.5)
        slow_play(FadeOut(density))

        # ------------------------------------------------------------------
        # Bar chart result — fully manual layout, no overlaps
        # ------------------------------------------------------------------
        chart_title = Text("What PrimeHelix measures", font_size=34)
        chart_title.to_edge(UP, buff=0.35)
        slow_play(Write(chart_title))

        # Move chart right to leave clean room for y-axis label.
        origin = LEFT * 3.2 + DOWN * 1.65
        chart_width = 7.4
        chart_height = 3.05

        x_axis = Line(origin, origin + RIGHT * chart_width, stroke_width=4)
        y_axis = Line(origin, origin + UP * chart_height, stroke_width=4)

        x_arrow = Triangle(fill_opacity=1).scale(0.12).rotate(-PI / 2)
        x_arrow.move_to(x_axis.get_end() + RIGHT * 0.12)

        y_arrow = Triangle(fill_opacity=1).scale(0.12)
        y_arrow.move_to(y_axis.get_end() + UP * 0.12)

        y_label = Text("Share of semiprimes", font_size=21).rotate(PI / 2)
        y_label.move_to(origin + LEFT * 1.0 + UP * (chart_height / 2))

        ticks = VGroup()
        tick_labels = VGroup()

        for pct in [20, 40, 60, 80]:
            y = origin[1] + (pct / 100) * chart_height
            tick = Line(
                start=[origin[0] - 0.11, y, 0],
                end=[origin[0] + 0.11, y, 0],
                stroke_width=3,
            )
            lab = Text(str(pct), font_size=18)
            lab.next_to(tick, LEFT, buff=0.18)
            ticks.add(tick)
            tick_labels.add(lab)

        slow_play(
            Create(x_axis),
            Create(y_axis),
            FadeIn(x_arrow),
            FadeIn(y_arrow),
            Write(y_label),
            FadeIn(ticks),
            FadeIn(tick_labels),
        )

        names = ["Lopsided", "Moderate", "Balanced"]
        values = [79, 20, 0.7]

        # Bar centers relative to chart origin.
        bar_positions = [1.0, 3.7, 6.3]
        bar_width = 0.75

        bars = VGroup()
        bar_labels = VGroup()
        value_labels = VGroup()

        for name, value, xpos in zip(names, values, bar_positions):
            height = max((value / 100) * chart_height, 0.035)
            center = origin + RIGHT * xpos + UP * (height / 2)

            bar = Rectangle(width=bar_width, height=height, stroke_width=4)
            bar.move_to(center)
            bars.add(bar)

            name_label = Text(name, font_size=21)
            name_label.move_to(origin + RIGHT * xpos + DOWN * 0.38)
            bar_labels.add(name_label)

            pct_text = "0.7%" if value < 1 else f"{value}%"
            pct_label = Text(pct_text, font_size=22)
            pct_label.next_to(bar, UP, buff=0.12)
            value_labels.add(pct_label)

        slow_play(*[GrowFromEdge(bar, DOWN) for bar in bars])
        slow_play(FadeIn(bar_labels), FadeIn(value_labels))
        slow_wait(1)

        # Fade the chart before showing the takeaway to avoid overlap.
        slow_play(
            FadeOut(chart_title),
            FadeOut(x_axis),
            FadeOut(y_axis),
            FadeOut(x_arrow),
            FadeOut(y_arrow),
            FadeOut(y_label),
            FadeOut(ticks),
            FadeOut(tick_labels),
            FadeOut(bars),
            FadeOut(bar_labels),
            FadeOut(value_labels),
        )

        takeaway = Text(
            "Balanced semiprimes are special — not typical.",
            font_size=34,
        )

        slow_play(Write(takeaway))
        slow_wait(1.5)
        slow_play(FadeOut(takeaway))

        # ------------------------------------------------------------------
        # Final message
        # ------------------------------------------------------------------
        final = Text(
            "Small primes divide often.\nLarge primes divide rarely.\nThat is why most semiprimes are lopsided.",
            font_size=36,
            line_spacing=0.9,
        )

        slow_play(Write(final))
        slow_wait(2)
