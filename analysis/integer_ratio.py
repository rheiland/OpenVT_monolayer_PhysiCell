def is_ratio_integer(a: float, b: float, tolerance: float = 1e-9) -> bool:
    """
    Test if the ratio a/b is an integer (within floating point tolerance).
    
    Args:
        a: numerator
        b: denominator
        tolerance: how close to an integer counts as integer (default 1e-9)
    
    Returns:
        True if a/b is (approximately) an integer, False otherwise
    """
    if b == 0:
        raise ValueError("Denominator cannot be zero")
    
    ratio = a / b
    return abs(ratio - round(ratio)) < tolerance


def check_ratio(a: float, b: float) -> None:
    try:
        ratio = a / b
        result = is_ratio_integer(a, b)
        status = "✓ INTEGER" if result else "✗ not integer"
        print(f"  {a} / {b} = {ratio:.10g}  →  {status}")
    except ValueError as e:
        print(f"  {a} / {b}  →  Error: {e}")


if __name__ == "__main__":
    print("=== Ratio Integer Tests ===\n")

    print("Exact integer ratios:")
    check_ratio(10.0, 2.0)       # 5.0
    check_ratio(7.5, 2.5)        # 3.0
    check_ratio(-9.0, 3.0)       # -3.0
    check_ratio(0.0, 5.0)        # 0.0

    print("\nNon-integer ratios:")
    check_ratio(1.0, 3.0)        # 0.333...
    check_ratio(7.0, 2.0)        # 3.5
    check_ratio(22.0, 7.0)       # ~3.1428...

    print("\nFloating point edge cases:")
    check_ratio(0.1 + 0.2, 0.3)  # Should be 1.0 (floating point noise)
    check_ratio(1e15, 1e5)        # 1e10
    check_ratio(1.0, 1e-10)       # 1e10

    print("\nError case:")
    check_ratio(5.0, 0.0)         # Division by zero
