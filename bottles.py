#!/usr/bin/env python3


def bottles(n):
    if n == 0:
        return "no more bottles"
    if n == 1:
        return "1 bottle"
    return f"{n} bottles"


def main():
    for n in range(99, 0, -1):
        print(f"{bottles(n)} of beer on the wall, {bottles(n)} of beer.")
        pronoun = "it" if n == 1 else "one"
        print(f"Take {pronoun} down and pass it around, {bottles(n - 1)} of beer on the wall.")
        print()
    print("No more bottles of beer on the wall, no more bottles of beer.")
    print("Go to the store and buy some more, 99 bottles of beer on the wall.")


if __name__ == "__main__":
    main()
