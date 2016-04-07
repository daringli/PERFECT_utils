#needed to place the offset text (the 10^x text by the axis) at arbitrary position: https://github.com/matplotlib/matplotlib/issues/4476

import types
import matplotlib.transforms as mtransforms
import numpy 
import matplotlib.pyplot as plt

def x_update_offset_text_position(self, bboxes, bboxes2):
    x, y = self.offsetText.get_position()

    if self.offset_text_position == 'bottom':
        if bboxes:
            bbox = mtransforms.Bbox.union(bboxes)
        else:
            bbox = self.axes.bbox
        y = bbox.ymin - self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        self.offsetText.set(va='top', ha='right')

    if self.offset_text_position == 'there':
        # y in axes coords, x i0n display coords
        self.offsetText.set_transform(mtransforms.blended_transform_factory(
                self.axes.transAxes, mtransforms.IdentityTransform()))

        top = self.axes.bbox.ymax
        right=self.axes.bbox.xmin
        #print x
        #print y
        y = top+15# - 2*self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        x = 1.03# + 0.1*self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    if self.offset_text_position == 'there2':
        # y in axes coords, x i0n display coords
        self.offsetText.set_transform(mtransforms.blended_transform_factory(
                self.axes.transAxes, mtransforms.IdentityTransform()))

        top = self.axes.bbox.ymax
        right=self.axes.bbox.xmin
        #print x
        #print y
        y = top+15# - 2*self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        x = 0.97# + 0.1*self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    else:
        if bboxes2:
            bbox = mtransforms.Bbox.union(bboxes2)
        else:
            bbox = self.axes.bbox
        y = bbox.ymax + self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        self.offsetText.set(va='bottom', ha='left')

    self.offsetText.set_position((x, y))

def y_update_offset_text_position(self, bboxes, bboxes2):
    x, y = self.offsetText.get_position()

    if self.offset_text_position == 'left':
        # y in axes coords, x in display coords
        self.offsetText.set_transform(mtransforms.blended_transform_factory(
                self.axes.transAxes, mtransforms.IdentityTransform()))

        top = self.axes.bbox.ymax
        y = top + self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    if self.offset_text_position == 'there':
        # y in axes coords, x in display coords
        self.offsetText.set_transform(mtransforms.blended_transform_factory(
                self.axes.transAxes, mtransforms.IdentityTransform()))

        top = self.axes.bbox.ymax
        right=self.axes.bbox.xmin
        #print x
        #print y
        y = top-10# - 2*self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        x = 0.03# + 0.1*self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    else:
        # x & y in display coords
        self.offsetText.set_transform(mtransforms.IdentityTransform())

        # Northwest of upper-right corner of right-hand extent of tick labels
        if bboxes2:
            bbox = mtransforms.Bbox.union(bboxes2)
        else:
            bbox = self.axes.bbox
        top, right = bbox.ymax, bbox.xmax
        x = right + self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        y = top + self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    self.offsetText.set_position((x, y))

def monkey_patch(axis,which):
    if which=='x':
        axis._update_offset_text_position = types.MethodType(x_update_offset_text_position, axis)
    if which=='y':
        axis._update_offset_text_position = types.MethodType(y_update_offset_text_position, axis)
   
