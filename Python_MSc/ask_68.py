# %%
import argparse

def closer(num,val):
    return min(val, key=lambda x:abs(x-num))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--parameter_1', type=int)
    parser.add_argument('--parameter_2',type = int, nargs='*')
    args = parser.parse_args()
    num_a = args.parameter_1
    val_a = args.parameter_2
    print(closer(num_a,val_a))
