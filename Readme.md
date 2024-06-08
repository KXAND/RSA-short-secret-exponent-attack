# Wiener Attack & Boneh-Durfee Attack

RSA 在私钥指数 $d$ 小于一定值的情况下是不安全的。

Wiener 攻击由 Wiener 在 "*Cryptanalysis of Short RSA Secret Exponents*" 提出，基于一种连分数分解算法。

Boneh-Durfee攻击由 Boneh 和 Durfee 在"*Cryptanalysis of RSA with private key d less than* $N^{0.292}$"提出，基于 Coppersmith 提出的一种格基归约算法（“格”是一种类似于向量空间的数学概念）。

Ciet 等人在 "*Short private exponent attacks on fast variants of RSA*"中讨论了多质数和 Takagi  Variant 中二者是否适用。

**这个仓库是对 Ciet 的论文的简单实现。具体来说就是：Takagi 的加解密、三质数和两质数 RSA 的 Wiener 攻击、Boneh-Durfee 攻击。**

- Wiener 攻击见[Multiprime](./Multiprime.py)、[CRT-RSA](./CRT-RSA.py)、[Takagi.py](./Takagi.py)。分别是经典两质数 RSA、多个（3个）质数 RSA 以及Takagi Variant 的情况。不过 Ciet 已经证明，Wiener 攻击是不适用于 Takagi Variant 的。这一部分是纯手搓的。

- Boneh-Durfee攻击见 [boneh_durfee](./boneh_durfee.py)。在最末调整注释情况调整为两质数或多质数的情况。这一部分是在 David Wong 的仓库[mimoo/RSA-and-LLL-attacks: attacking RSA via lattice reductions (LLL)](https://github.com/mimoo/RSA-and-LLL-attacks)上修改而来。同时非常推荐观看David Wong 提供的讲解视频。

## 代码环境

代码使用了 SageMath 提供的数学库，所以无法在 Windows 系统直接运行。可以使用 VScode 的远程连接加 WSL 使用。

运行

```bash
apt install sagemath
```
安装 sagemath。（也许可以只运行`apt install python3-sage`，但是我不确定）

我使用的 Python 版本：3.10

我使用的Sage 版本：9.5
