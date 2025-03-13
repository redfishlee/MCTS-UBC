import math
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
from MEACRNG.MolEncoder.C1MolLib import C1MolLib as c
from PIL import Image, ImageFilter, ImageDraw, ImageFont# pillow库，python3中一个处理图像的工具
from PIL import ImageOps
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

'''
rdkit不支持绘制过渡态的键，现在将芳香键改成过渡态键的图形，不再支持过渡态+芳香键，因芳香键很少出现
修改位置：
这段代码修改自rdkit/Chem.Draw.MolDrawing
将_drawBond中对于芳香键的部分进行了修改

elif bType == Chem.BondType.AROMATIC:
     #addDefaultLine(newpos, newnbrPos)
     fp1, fp2 = self._offsetDblBond(newpos, newnbrPos, bond, atom, nbr, conf,direction=0)
     addDefaultLine(fp1, fp2, dash=self.drawingOptions.dash)

原来为：

elif bType == Chem.BondType.AROMATIC:
     addDefaultLine(newpos, newnbrPos)
     fp1, fp2 = self._offsetDblBond(newpos, newnbrPos, bond, atom, nbr, conf)
     addDefaultLine(fp1, fp2, dash=self.drawingOptions.dash)

因为是采用芳香键来代替过渡态键，因此修改
可以将本文件夹下的 ./MolDrawing_mod_from_rdkit重命名为MolDrawing，替换rdkit的同名文件


'''


def __turn_to_name_with_sub(name):# 使用LaTeX公式修正数字下标，用于matplotlib作图
    new_name = ""
    for i in name:
        assert isinstance(i, str)
        if i.isdigit():# 检测字符串i是否完全由数字组成
            new_name += "$_{%s}$" % i# $_{}$是LaTeX的公式，表示下标
        else:
            new_name += i
    return new_name


def show_mol(mol, regenerate=True, atom_fontsize=12, line_width=1.2, **kwargs):# 绘制分子图模型
    if regenerate:
        mol = Chem.MolFromSmarts(Chem.MolToSmarts(mol))# 印象里是一个编码问题，需要再生一次才能正常显示
    # 修改绘图设置.这些都是绘图的参数，具体可查询rdkit的官方文档
    Draw.DrawingOptions.atomLabelFontSize = atom_fontsize# 原子符号大小
    Draw.DrawingOptions.bondLineWidth = line_width# 键线宽度
    Draw.DrawingOptions.dash = (8, 8)# 虚线长度

    Draw.DrawingOptions.elemDict[1] = (0,0,1)# 调整某些原子颜色.原子序数:(红，绿，蓝).1是H，0不知道是什么

    img = Draw.MolToImage(mol, **kwargs)# 可以多个mol一起逐个绘图
    img.show()


def show_mol_array(mol_list, species_name_replace=None):# 绘制图阵
    Draw.DrawingOptions.atomLabelFontSize = 27
    Draw.DrawingOptions.bondLineWidth = 4
    Draw.DrawingOptions.dash = (8, 8)

    print("注意把修改的内容提交到原来的MEACRNG repo中")
    mol_num = len(mol_list)
    # 开方计算行列
    n_col = int(math.sqrt(mol_num))
    if n_col * n_col < mol_num:
        fig, axes = plt.subplots(n_col + 1, n_col, figsize=(9, 9))
    else:
        fig, axes = plt.subplots(n_col, n_col, figsize=(9, 9))

    plt.subplots_adjust(wspace=0, hspace=0)

    axes = axes.reshape(-1)

    for i in range(axes.shape[0]):
        axes[i].axis('off')

    for i in range(len(axes)):

        # imshow 图像
        # 得到图像
        img = Draw.MolToImage(mol_list[i], size=(600, 600), fitImage=True)
        # 检测白边，裁去白边
        ivt_image = ImageOps.invert(img)
        bbox = ivt_image.getbbox()
        img = img.crop([bbox[0], bbox[1], bbox[2], bbox[3]])
        # 然后把一边放大到300
        max_ = max(img.width, img.height)
        ratio = 300 / max_
        img = img.resize((int(img.width * ratio), int(img.height * ratio)))
        Image._show(img)

        # img = img.filter(ImageFilter.BLUR)
        # 加上白色背景，使其fix大小
        bg = Image.new("RGB", (310, 310), color="#FFFFFF")
        bw, bh = bg.size
        lw, lh = img.size
        # paste的起点，用于计算图例文字的起点
        start_p = (int((bw - lw) / 2), int((bh - lh) / 2))
        bg.paste(img, start_p)
        axes[i].imshow(bg)

        label_font_size = 15
        if species_name_replace is not None:
            axes[i].text(0, 0, __turn_to_name_with_sub(species_name_replace[i]), fontsize=label_font_size)

        else:
            axes[i].set_title(Chem.MolToSmiles(mol_list[i]), fontsize=label_font_size)

    fig.tight_layout()  # 调整整体空白
    plt.subplots_adjust(wspace=0, hspace=.1)  # 调整子图间距

    plt.show()


def show_mol_array_using_PIL(mol_list, species_name_replace=None,save_filename="fig.png"):
    Draw.DrawingOptions.atomLabelFontSize = 27
    Draw.DrawingOptions.bondLineWidth = 4
    Draw.DrawingOptions.dash = (8, 8)
    #Draw.DrawingOptions.elemDict[1] = (0,0,1)

    print("注意把修改的内容提交到原来的MEACRNG repo中")
    mol_num = len(mol_list)
    # 开方计算行列
    n_col = int(math.sqrt(mol_num))
    if n_col * n_col < mol_num:
        x, y = n_col + 1, n_col
    else:
        x, y = n_col, n_col

    now_fig_index = -1
    max_ = max(x, y)
    fontStyle = ImageFont.truetype("C:/Windows/Fonts/ARIAL.ttf", size=30)
    fontStyle_for_sub = ImageFont.truetype("C:/Windows/Fonts/ARIAL.ttf", size=20)
    total_bg = Image.new("RGB", (320 * max_, 320 * max_), color="#FFFFFF")
    draw = ImageDraw.Draw(total_bg)
    for _x in range(x):
        for _y in range(y):
            now_fig_index += 1
            if now_fig_index > mol_num - 1: continue
            # imshow 图像
            # 得到图像
            img = Draw.MolToImage(mol_list[now_fig_index], size=(600, 600), fitImage=True)
            # 检测白边，裁去白边
            ivt_image = ImageOps.invert(img)
            bbox = ivt_image.getbbox()
            img = img.crop([bbox[0]-10, bbox[1]-10, bbox[2] + 10, bbox[3] + 10])

            # img = img.filter(ImageFilter.BLUR)
            # 加上白色背景，使其fix大小
            bg = Image.new("RGB", (320, 320), color="#FFFFFF")
            bw, bh = bg.size
            lw, lh = img.size
            # paste的起点，用于计算图例文字的起点
            start_p = (int((bw - lw) / 2), int((bh - lh) / 2))
            bg.paste(img, start_p)
            # 把这个图片贴到大bg上，要求居中
            sp_in_total_bg = (_x * 310+10, _y * 310+10)
            total_bg.paste(bg, sp_in_total_bg)
            # 带下标的draw方法
            now_char_pos = 0
            for i in range(len(species_name_replace[now_fig_index])):
                char = species_name_replace[now_fig_index][i]
                assert isinstance(char, str)
                if char.isdigit():
                    draw.text(xy=(sp_in_total_bg[0] + now_char_pos, sp_in_total_bg[1]+15),
                              text=char,
                              fill=(0, 0, 0), font=fontStyle_for_sub)
                    now_char_pos += fontStyle_for_sub.size * 0.75
                elif char in ["(",")","-"]:
                    draw.text(xy=(sp_in_total_bg[0] + now_char_pos, sp_in_total_bg[1]),
                              text=char,
                              fill=(0, 0, 0), font=fontStyle)
                    now_char_pos += fontStyle.size * 0.4

                else:
                    draw.text(xy=(sp_in_total_bg[0] + now_char_pos, sp_in_total_bg[1]),
                              text=char,
                              fill=(0, 0, 0), font=fontStyle)
                    now_char_pos += fontStyle.size * 0.75

            # label_font_size = 15
            # if species_name_replace is not None:
            #     axes[i].text(0, 0, __turn_to_name_with_sub(species_name_replace[i]), fontsize=label_font_size)
            #
            # else:
            #     axes[i].set_title(Chem.MolToSmiles(mol_list[i]), fontsize=label_font_size)
    total_bg.save(save_filename)
    Image._show(total_bg)


def tst_show_mol_array():
    show_mol_array([
        c.CHO,
        c.CH3,
        c.CH,
        c.COOH_single_CO_bond,
        c.COOH_single_CO_bond,
        c.COOH_single_CO_bond,
        c.COOH_single_CO_bond,
        c.COOH_single_CO_bond,
    ])


if __name__ == '__main__':
    tst_show_mol_array()
