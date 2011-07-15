object MainWindow: TMainWindow
  Left = 338
  Top = 219
  Width = 819
  Height = 550
  Caption = #1048#1085#1090#1077#1075#1088#1080#1088#1086#1074#1072#1085#1080#1077' '#1082#1072#1089#1082#1072#1076#1085#1086'-'#1089#1077#1090#1086#1095#1085#1099#1084' '#1084#1077#1090#1086#1076#1086#1084
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  OnDestroy = FormDestroy
  DesignSize = (
    811
    516)
  PixelsPerInch = 96
  TextHeight = 13
  object Notebook: TTabbedNotebook
    Left = 0
    Top = -24
    Width = 981
    Height = 532
    Anchors = [akLeft, akTop, akRight, akBottom]
    TabFont.Charset = DEFAULT_CHARSET
    TabFont.Color = clBtnText
    TabFont.Height = -11
    TabFont.Name = 'MS Sans Serif'
    TabFont.Style = []
    TabOrder = 0
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Initiation'
      object GroupOutput: TGroupBox
        Left = 296
        Top = 24
        Width = 337
        Height = 313
        Caption = #1048#1084#1077#1085#1072' '#1092#1072#1081#1083#1086#1074' '#1074#1099#1074#1086#1076#1072
        TabOrder = 0
        object Label1: TLabel
          Left = 16
          Top = 36
          Width = 89
          Height = 13
          Caption = #1060#1072#1081#1083' '#1089#1086#1086#1073#1097#1077#1085#1080#1081
        end
        object Label2: TLabel
          Left = 16
          Top = 64
          Width = 127
          Height = 13
          Caption = #1042#1099#1074#1086#1076' '#1079#1085#1072#1095#1077#1085#1080#1081' '#1101#1085#1077#1088#1075#1080#1080
        end
        object Label3: TLabel
          Left = 16
          Top = 96
          Width = 133
          Height = 13
          Caption = #1042#1099#1074#1086#1076' '#1079#1085#1072#1095#1077#1085#1080#1081' '#1089#1082#1086#1088#1086#1089#1090#1080
        end
        object Label4: TLabel
          Left = 16
          Top = 156
          Width = 157
          Height = 13
          Caption = #1042#1099#1074#1086#1076' '#1082#1072#1089#1082#1072#1076#1085#1099#1093' '#1087#1077#1088#1077#1084#1077#1085#1085#1099#1093
        end
        object Label5: TLabel
          Left = 16
          Top = 196
          Width = 70
          Height = 13
          Caption = #1060#1072#1081#1083' '#1086#1096#1080#1073#1086#1082
        end
        object Label6: TLabel
          Left = 16
          Top = 260
          Width = 64
          Height = 13
          Caption = #1060#1072#1081#1083' '#1076#1072#1084#1087#1072
        end
        object Label7: TLabel
          Left = 16
          Top = 124
          Width = 133
          Height = 13
          Caption = #1042#1099#1074#1086#1076' '#1079#1085#1072#1095#1077#1085#1080#1081' '#1074#1103#1079#1082#1086#1089#1090#1080
        end
        object EditMessage: TEdit
          Left = 152
          Top = 32
          Width = 145
          Height = 21
          TabOrder = 0
          Text = 'message.dat'
          OnKeyUp = EditMessageKeyUp
        end
        object EditEnergy: TEdit
          Left = 152
          Top = 60
          Width = 145
          Height = 21
          TabOrder = 1
          Text = 'energy.dat'
          OnKeyUp = EditEnergyKeyUp
        end
        object EditNu: TEdit
          Left = 152
          Top = 120
          Width = 145
          Height = 21
          TabOrder = 3
          Text = 'nut.dat'
          OnKeyUp = EditNuKeyUp
        end
        object EditKaskvar: TEdit
          Left = 184
          Top = 152
          Width = 113
          Height = 21
          TabOrder = 4
          Text = 'kaskvar.dat'
          OnKeyUp = EditKaskvarKeyUp
        end
        object EditError: TEdit
          Left = 152
          Top = 192
          Width = 145
          Height = 21
          TabOrder = 5
          Text = 'error.err'
          OnKeyUp = EditErrorKeyUp
        end
        object EditDump: TEdit
          Left = 152
          Top = 256
          Width = 145
          Height = 21
          TabOrder = 6
          Text = 'dump.dat'
          OnChange = EditDumpChange
        end
        object EditV: TEdit
          Left = 152
          Top = 92
          Width = 145
          Height = 21
          TabOrder = 2
          Text = 'vv.dat'
          OnKeyUp = EditVKeyUp
        end
      end
      object GroupInit: TGroupBox
        Left = 24
        Top = 24
        Width = 249
        Height = 233
        Caption = #1042#1074#1086#1076' '#1085#1072#1095#1072#1083#1100#1085#1099#1093' '#1076#1072#1085#1085#1099#1093
        TabOrder = 1
        object Label8: TLabel
          Left = 32
          Top = 40
          Width = 60
          Height = 13
          Caption = #1048#1084#1103' '#1092#1072#1081#1083#1072':'
        end
        object EditInit: TEdit
          Left = 32
          Top = 64
          Width = 97
          Height = 21
          TabOrder = 0
          Text = 'init.dat'
          OnKeyUp = EditInitKeyUp
        end
        object InputInit: TButton
          Left = 32
          Top = 168
          Width = 169
          Height = 41
          Caption = #1042#1074#1077#1089#1090#1080' '#1076#1072#1085#1085#1099#1077' '#1074' '#1087#1088#1086#1075#1088#1072#1084#1084#1091
          TabOrder = 1
          OnClick = InputInitClick
        end
        object RadioInitVel: TRadioButton
          Left = 32
          Top = 104
          Width = 169
          Height = 17
          Caption = #1042#1074#1077#1089#1090#1080' '#1085#1072#1095#1072#1083#1100#1085#1091#1102' '#1089#1082#1086#1088#1086#1089#1090#1100
          Checked = True
          TabOrder = 2
          TabStop = True
        end
        object RadioInitContinue: TRadioButton
          Left = 32
          Top = 136
          Width = 113
          Height = 17
          Caption = #1055#1088#1086#1076#1086#1083#1078#1080#1090#1100' '#1089#1095#1077#1090
          TabOrder = 3
        end
      end
    end
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Calculation'
      object BeginWork: TButton
        Left = 480
        Top = 56
        Width = 145
        Height = 45
        Caption = #1053#1072#1095#1072#1090#1100' '#1089#1095#1077#1090
        ParentShowHint = False
        ShowHint = False
        TabOrder = 1
        OnClick = BeginWorkClick
      end
      object Screen: TMemo
        Left = 24
        Top = 48
        Width = 417
        Height = 329
        Color = clCaptionText
        ScrollBars = ssVertical
        TabOrder = 0
      end
      object DoDump: TButton
        Left = 480
        Top = 227
        Width = 145
        Height = 45
        Caption = #1057#1073#1088#1086#1089#1080#1090#1100' '#1076#1072#1085#1085#1099#1077' '#1074' '#1092#1072#1081#1083
        TabOrder = 2
        OnClick = DoDumpClick
      end
      object Exit: TButton
        Left = 480
        Top = 312
        Width = 145
        Height = 45
        Caption = #1054#1082#1086#1085#1095#1072#1085#1080#1077' '#1088#1072#1073#1086#1090#1099
        TabOrder = 3
        OnClick = ExitClick
      end
      object ButtonRecontinue: TButton
        Left = 480
        Top = 141
        Width = 145
        Height = 45
        Caption = #1053#1072#1095#1072#1090#1100' '#1089#1095#1077#1090' '#1089#1099#1079#1085#1086#1074#1072
        TabOrder = 4
        Visible = False
        OnClick = RecontinueClick
      end
      object ClearScr: TButton
        Left = 24
        Top = 392
        Width = 65
        Height = 25
        Caption = #1054#1095#1080#1089#1090#1080#1090#1100
        TabOrder = 5
        OnClick = ClearScrClick
      end
      object SelPrior: TComboBox
        Left = 656
        Top = 152
        Width = 113
        Height = 21
        Hint = #1087#1086' '#1091#1084#1086#1083#1095#1072#1085#1080#1102' '#1089#1072#1084#1099#1081' '#1085#1080#1079#1082#1080#1081
        ItemHeight = 13
        ParentShowHint = False
        ShowHint = True
        TabOrder = 6
        Text = #1055#1088#1080#1086#1088#1080#1090#1077#1090
        OnChange = SelPriorChange
        Items.Strings = (
          '0.'#1074' '#1073#1077#1079#1076#1077#1083#1100#1077
          '1.'#1085#1080#1079#1082#1080#1081
          '2.'#1085#1080#1078#1077' '#1089#1088#1077#1076#1085#1077#1075#1086
          '3.'#1085#1086#1088#1084#1072#1083#1100#1085#1099#1081
          '4.'#1074#1099#1096#1077' '#1089#1088#1077#1076#1085#1077#1075#1086
          '5.'#1074#1099#1089#1086#1082#1080#1081
          '6.'#1088#1077#1072#1083#1100#1085#1086#1075#1086' '#1074#1088#1077#1084#1077#1085#1080)
      end
    end
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Printing'
      object Label9: TLabel
        Left = 384
        Top = 56
        Width = 184
        Height = 13
        Caption = #1050#1072#1082' '#1095#1072#1089#1090#1086' '#1087#1088#1086#1074#1077#1088#1103#1090#1100' '#1088#1072#1089#1093#1086#1076#1080#1084#1086#1089#1090#1100':'
      end
      object Label10: TLabel
        Left = 416
        Top = 144
        Width = 40
        Height = 13
        Caption = #1082#1072#1078#1076#1099#1081
      end
      object Label11: TLabel
        Left = 624
        Top = 85
        Width = 81
        Height = 13
        Caption = #1096#1072#1075' '#1087#1086' '#1074#1088#1077#1084#1077#1085#1080
      end
      object Label12: TLabel
        Left = 624
        Top = 144
        Width = 81
        Height = 13
        Caption = #1096#1072#1075' '#1087#1086' '#1074#1088#1077#1084#1077#1085#1080
      end
      object Label13: TLabel
        Left = 416
        Top = 85
        Width = 40
        Height = 13
        Caption = #1082#1072#1078#1076#1099#1081
      end
      object Label14: TLabel
        Left = 384
        Top = 120
        Width = 208
        Height = 13
        Caption = #1050#1072#1082' '#1095#1072#1089#1090#1086' '#1074#1099#1074#1086#1076#1080#1090#1100' '#1085#1072' '#1101#1082#1088#1072#1085' '#1080' '#1074' '#1092#1072#1081#1083#1099':'
      end
      object CheckScreenOutput: TCheckBox
        Left = 72
        Top = 88
        Width = 169
        Height = 17
        Hint = #1091#1076#1086#1073#1085#1086
        Caption = #1042#1099#1074#1086#1076#1080#1090#1100' '#1076#1072#1085#1085#1099#1077' '#1085#1072' '#1101#1082#1088#1072#1085
        Checked = True
        ParentShowHint = False
        ShowHint = True
        State = cbChecked
        TabOrder = 1
      end
      object CheckFileOutput: TCheckBox
        Left = 72
        Top = 56
        Width = 177
        Height = 17
        Hint = #1073#1077#1079' '#1101#1090#1086#1075#1086' '#1085#1077#1083#1100#1079#1103
        Caption = #1042#1099#1074#1086#1076#1080#1090#1100' '#1076#1072#1085#1085#1099#1077' '#1074' '#1092#1072#1081#1083#1099
        Checked = True
        ParentShowHint = False
        ShowHint = True
        State = cbChecked
        TabOrder = 0
      end
      object CheckClear: TCheckBox
        Left = 72
        Top = 120
        Width = 249
        Height = 17
        Caption = #1054#1095#1080#1097#1072#1090#1100' '#1101#1082#1088#1072#1085' '#1087#1086#1089#1083#1077' '#1074#1099#1074#1086#1076#1072' '#1086#1076#1085#1086#1075#1086' '#1096#1072#1075#1072
        Checked = True
        State = cbChecked
        TabOrder = 2
      end
      object EditCheckStep: TEdit
        Left = 480
        Top = 80
        Width = 120
        Height = 21
        TabOrder = 4
        Text = '10'
        OnChange = EditCheckStepChanging
      end
      object UpDownCheck: TUpDown
        Left = 600
        Top = 80
        Width = 15
        Height = 21
        Associate = EditCheckStep
        Min = 1
        Position = 10
        TabOrder = 3
        Wrap = False
      end
      object EditOutStep: TEdit
        Left = 480
        Top = 144
        Width = 120
        Height = 21
        TabOrder = 5
        Text = '20'
        OnChange = EditOutStepChanging
      end
      object UpDownOut: TUpDown
        Left = 600
        Top = 144
        Width = 15
        Height = 21
        Associate = EditOutStep
        Min = 1
        Max = 1000
        Position = 20
        TabOrder = 6
        Wrap = False
      end
    end
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Parameters'
      object Label15: TLabel
        Left = 96
        Top = 68
        Width = 20
        Height = 13
        Caption = 'Re='
      end
      object GroupTurbulent: TGroupBox
        Left = 40
        Top = 264
        Width = 425
        Height = 153
        Caption = #1058#1091#1088#1073#1091#1083#1077#1085#1090#1085#1072#1103' '#1074#1103#1079#1082#1086#1089#1090#1100
        TabOrder = 0
        object Label18: TLabel
          Left = 24
          Top = 60
          Width = 180
          Height = 13
          Caption = #1050#1086#1101#1092#1092#1080#1094#1080#1077#1085#1090' '#1087#1088#1086#1087#1086#1088#1094#1080#1086#1085#1072#1083#1100#1085#1086#1089#1090#1080
        end
        object CheckNut: TCheckBox
          Left = 16
          Top = 32
          Width = 369
          Height = 17
          Caption = #1048#1089#1087#1086#1083#1100#1079#1086#1074#1072#1090#1100' '#1090#1091#1088#1073#1091#1083#1077#1085#1090#1085#1091#1102' '#1074#1103#1079#1082#1086#1089#1090#1100' '#1095#1077#1088#1077#1079' '#1082#1072#1089#1082#1072#1076#1085#1099#1077' '#1092#1091#1085#1082#1094#1080#1080
          Checked = True
          State = cbChecked
          TabOrder = 0
          OnClick = CheckNutClick
        end
        object EditMaschtab: TEdit
          Left = 224
          Top = 56
          Width = 121
          Height = 21
          TabOrder = 1
          Text = '1000'
          OnChange = EditMaschtabChange
        end
      end
      object EditRe: TEdit
        Left = 120
        Top = 64
        Width = 121
        Height = 21
        TabOrder = 1
        Text = '10000'
        OnKeyUp = EditReKeyUp
      end
      object GroupInitVariant: TGroupBox
        Left = 40
        Top = 104
        Width = 305
        Height = 145
        Caption = #1042#1099#1073#1077#1088#1080#1090#1077' '#1090#1080#1087' '#1085#1072#1095#1072#1083#1100#1085#1099#1093' '#1091#1089#1083#1086#1074#1080#1081' '#1076#1083#1103' '#1089#1082#1086#1088#1086#1089#1090#1080
        TabOrder = 2
        object RadioZeroInit: TRadioButton
          Left = 24
          Top = 24
          Width = 113
          Height = 17
          Caption = #1053#1091#1083#1077#1074#1072#1103
          TabOrder = 0
        end
        object RadioParabole: TRadioButton
          Left = 24
          Top = 56
          Width = 129
          Height = 17
          Caption = #1055#1088#1086#1092#1080#1083#1100' '#1055#1091#1072#1079#1077#1081#1083#1103
          Checked = True
          TabOrder = 1
          TabStop = True
        end
      end
      object GroupNoise: TGroupBox
        Left = 496
        Top = 264
        Width = 217
        Height = 153
        Caption = #1055#1072#1088#1072#1084#1077#1090#1088#1099' '#1096#1091#1084#1072
        TabOrder = 3
        object Label16: TLabel
          Left = 24
          Top = 56
          Width = 78
          Height = 13
          Caption = #1042#1077#1083#1080#1095#1080#1085#1072' '#1096#1091#1084#1072
        end
        object Label17: TLabel
          Left = 24
          Top = 116
          Width = 80
          Height = 13
          Caption = #1042#1077#1083#1080#1095#1080#1085#1072' '#1074#1086#1079#1084'.'
        end
        object CheckNoise: TCheckBox
          Left = 24
          Top = 24
          Width = 97
          Height = 17
          Caption = #1042#1082#1083#1102#1095#1080#1090#1100' '#1096#1091#1084
          Checked = True
          State = cbChecked
          TabOrder = 0
          OnClick = CheckNoiseClick
        end
        object EditNoise: TEdit
          Left = 120
          Top = 52
          Width = 65
          Height = 21
          TabOrder = 1
          Text = '0.1'
          OnKeyUp = EditNoiseKeyUp
        end
        object CheckNorm: TCheckBox
          Left = 24
          Top = 88
          Width = 177
          Height = 17
          Caption = #1042#1082#1083#1102#1095#1080#1090#1100' '#1085#1086#1088#1084'. '#1074#1086#1079#1084#1091#1097#1077#1085#1080#1103
          TabOrder = 2
          OnClick = CheckNormClick
        end
        object EditNorm: TEdit
          Left = 120
          Top = 112
          Width = 65
          Height = 21
          Enabled = False
          TabOrder = 3
          Text = '0'
          OnChange = EditNormChange
        end
      end
      object GroupGrid: TGroupBox
        Left = 392
        Top = 104
        Width = 353
        Height = 145
        Caption = #1055#1072#1088#1072#1084#1077#1090#1088#1099' '#1089#1077#1090#1082#1080
        TabOrder = 4
        object Label19: TLabel
          Left = 16
          Top = 24
          Width = 89
          Height = 13
          Caption = #1050#1086#1083'-'#1074#1086' '#1091#1079#1083#1086#1074' '#1087#1086' x'
        end
        object Label20: TLabel
          Left = 16
          Top = 100
          Width = 89
          Height = 13
          Caption = #1050#1086#1083'-'#1074#1086' '#1091#1079#1083#1086#1074' '#1087#1086' z'
        end
        object Label21: TLabel
          Left = 16
          Top = 62
          Width = 89
          Height = 13
          Caption = #1050#1086#1083'-'#1074#1086' '#1091#1079#1083#1086#1074' '#1087#1086' y'
        end
        object EditNx: TEdit
          Left = 144
          Top = 20
          Width = 65
          Height = 21
          TabOrder = 0
          Text = 'EditNx'
        end
        object EditNy: TEdit
          Left = 144
          Top = 59
          Width = 65
          Height = 21
          TabOrder = 1
          Text = 'EditNy'
        end
        object EditNz: TEdit
          Left = 144
          Top = 96
          Width = 65
          Height = 21
          TabOrder = 2
          Text = 'EditNz'
        end
        object ButtonSaveNi: TButton
          Left = 224
          Top = 45
          Width = 113
          Height = 49
          Caption = #1055#1088#1080#1085#1103#1090#1100' '#1080#1079#1084#1077#1085#1077#1085#1080#1103
          TabOrder = 3
          OnClick = ButtonSaveNiClick
        end
      end
    end
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Tasks'
    end
  end
  object StatusBarCalc: TStatusBar
    Left = 0
    Top = 496
    Width = 811
    Height = 20
    Panels = <
      item
        Text = #1053#1072#1095#1072#1083#1086' '#1088#1072#1073#1086#1090#1099' ('#1085#1072#1078#1084#1080#1090#1077' '#1082#1085#1086#1087#1082#1091')'
        Width = 200
      end
      item
        Text = #1053#1080#1095#1077#1075#1086' '#1085#1077' '#1089#1076#1077#1083#1072#1085#1086
        Width = 50
      end>
    SimplePanel = False
  end
end
