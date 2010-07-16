object Form1: TForm1
  Left = 215
  Top = 107
  Width = 696
  Height = 480
  Caption = 'Form1'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Button1: TButton
    Left = 216
    Top = 80
    Width = 75
    Height = 25
    Caption = 'Button1'
    TabOrder = 0
    OnClick = Button1Click
  end
  object Button2: TButton
    Left = 336
    Top = 80
    Width = 75
    Height = 25
    Caption = 'Button2'
    TabOrder = 1
    OnKeyPress = Button2KeyPress
  end
  object ApplicationEvents1: TApplicationEvents
    Left = 96
    Top = 40
  end
end
